"""
agent_orchestrator.py — the multi-agent coordinator.

Flow for each question:
  1. route()        — pick the 1-3 relevant specialists (agents.py)
  2. run specialists in parallel threads — each sees the SAME grounded context
     (retrieval.get_protein_context) + dataset-tool results when relevant
  3. synthesize     — merge specialist answers into one cited answer + an
     "Agents consulted" trace
  4. stream the final answer

Reuses chat_backend's key/client resolution and prompts.format_context_prompt so
grounding and citations are identical to the single-protein chat. Offline (no API
key) falls back to a structured multi-specialist report card.
"""
from __future__ import annotations
import os
import concurrent.futures as cf
from typing import Any, Generator, Optional

try:
    import agents
    import analysis_tools
    from prompts import format_context_prompt
    from chat_backend import make_client, resolve_api_key, format_report_card, DEFAULT_MODEL
except ImportError:  # package-relative
    from . import agents, analysis_tools  # type: ignore
    from .prompts import format_context_prompt  # type: ignore
    from .chat_backend import make_client, resolve_api_key, format_report_card, DEFAULT_MODEL  # type: ignore

SPECIALIST_MAX_TOKENS = 900
SYNTHESIS_MAX_TOKENS = 1400


# ---------------------------------------------------------------------------
# Dataset tool gathering — the Data Analyst gets precomputed cross-cutting data.
# ---------------------------------------------------------------------------
def _gather_dataset_results(question: str, uid: Optional[str]) -> dict[str, Any]:
    """Run the cheap dataset tools relevant to the question. Returns
    {tool_name: result_dict}. All are precomputed-table lookups (fast)."""
    out: dict[str, Any] = {}
    q = question.lower()
    try:
        if uid:
            out["compare_conditions"] = analysis_tools.compare_conditions(uid)
        if agents.is_dataset_level(question):
            if any(k in q for k in ["aging-specific", "aging specific", "specific", "vs stress", "generic"]):
                out["aging_specific_hits"] = analysis_tools.aging_specific_hits(n=15)
            if any(k in q for k in ["top", "most", "largest", "strongest", "mover"]):
                # default to day9 unless another condition named
                cond = "day9"
                for c in analysis_tools.COND_LABEL:
                    if c in q or analysis_tools.COND_LABEL[c].lower() in q:
                        cond = c; break
                out["top_movers"] = analysis_tools.top_movers(cond, n=15)
            if any(k in q for k in ["enrich", "buried", "structur", "near site", "functional site"]):
                out["structural_enrichment"] = analysis_tools.structural_enrichment_summary()
            if any(k in q for k in ["disease", "candidate", "discovery", "variant"]):
                out["discovery_top"] = analysis_tools.discovery_top(n=10)
    except Exception as e:
        out["_error"] = f"{type(e).__name__}: {e}"
    return out


def _dataset_block(results: dict[str, Any]) -> str:
    """Render dataset-tool results as a text block appended to the Data Analyst's context."""
    if not results:
        return ""
    lines = ["\nDATASET ANALYSIS (precomputed, whole-dataset — cite as [S5]/[S6]):"]
    for name, r in results.items():
        if name.startswith("_") or not isinstance(r, dict):
            continue
        if r.get("summary"):
            lines.append(f"  • {name}: {r['summary']}")
    return "\n".join(lines) if len(lines) > 1 else ""


# ---------------------------------------------------------------------------
# Specialist call
# ---------------------------------------------------------------------------
def _run_specialist(agent_id: str, question: str, ctx: dict[str, Any],
                    dataset_block: str, client, model: str) -> dict[str, Any]:
    """Call one specialist. Returns {id, name, icon, text, error}."""
    spec = agents.AGENTS[agent_id]
    user_msg = format_context_prompt(ctx, question)
    if dataset_block and agent_id == "data_analyst":
        user_msg += "\n" + dataset_block
    try:
        resp = client.messages.create(
            model=model,
            system=spec["system_prompt"],
            max_tokens=SPECIALIST_MAX_TOKENS,
            messages=[{"role": "user", "content": user_msg}],
        )
        text = "".join(b.text for b in resp.content if getattr(b, "type", None) == "text")
        return {"id": agent_id, "name": spec["name"], "icon": spec["icon"], "text": text.strip(),
                "usage": {"in": resp.usage.input_tokens, "out": resp.usage.output_tokens}}
    except Exception as e:
        return {"id": agent_id, "name": spec["name"], "icon": spec["icon"],
                "text": "", "error": f"{type(e).__name__}: {e}"}


_SYNTHESIS_SYSTEM = """You are the coordinator of a panel of specialist agents analyzing a C. elegans LiP-MS aging dataset. You are given the user's question and each specialist's grounded, cited response. Produce ONE integrated answer.

Rules:
- Synthesize — do not just concatenate. Resolve overlaps, order the reasoning logically, and lead with the direct answer to the question.
- PRESERVE every bracketed citation [S1]-[S8] from the specialists; do not invent new ones or drop them. Keep the evidence-type discipline (measured [S5] vs ML [S6] vs homology [S7] vs AlphaFold-predicted [S8]).
- If specialists disagree or add different angles, integrate both. If one specialist found no relevant data, don't pad — just omit.
- Be concise and well-structured: a direct answer, then supporting detail. 3-6 short paragraphs or grouped bullets.
- Do not mention "the specialists" or the panel machinery in the answer body — just give the integrated scientific answer.
"""


def _synthesize_stream(question: str, specialist_outputs: list[dict], ctx: dict,
                       client, model: str) -> Generator[str, None, dict]:
    """Stream the synthesized final answer merging specialist outputs."""
    contributions = []
    for o in specialist_outputs:
        if o.get("text"):
            contributions.append(f"### {o['name']} ({', '.join(agents.AGENTS[o['id']]['sources'])}):\n{o['text']}")
    joined = "\n\n".join(contributions) if contributions else "(no specialist produced a grounded answer)"
    ident = ctx.get("identity", {})
    prot = f"{ident.get('gene_symbol','?')} ({ident.get('uniprot_id','?')})" if ident else "the dataset"
    synth_user = (f"USER QUESTION: {question}\n\nProtein in focus: {prot}\n\n"
                  f"SPECIALIST RESPONSES TO INTEGRATE:\n\n{joined}")
    try:
        with client.messages.stream(
            model=model, system=_SYNTHESIS_SYSTEM,
            max_tokens=SYNTHESIS_MAX_TOKENS,
            messages=[{"role": "user", "content": synth_user}],
        ) as stream:
            for t in stream.text_stream:
                yield t
            final = stream.get_final_message()
            yield {"usage": {"input_tokens": final.usage.input_tokens,
                             "output_tokens": final.usage.output_tokens},
                   "model": final.model, "stop_reason": final.stop_reason}
    except Exception as e:
        yield f"\n\n**Synthesis error:** `{type(e).__name__}: {e}`"
        yield {"usage": None, "model": model, "error": str(e)}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def ask_agents(question: str, ctx: dict[str, Any],
               uid: Optional[str] = None,
               model: str = DEFAULT_MODEL,
               api_key: Optional[str] = None,
               max_specialists: int = 3) -> Generator[Any, None, dict]:
    """Route → run specialists in parallel → synthesize. Streaming generator.

    Yields:
      - {"event": "route", "specialists": [...], "scores": {...}}   (once, first)
      - {"event": "specialist_done", "id", "name", "icon", "ok"}    (per specialist)
      - str chunks of the synthesized answer
      - final dict {"event":"final", "trace":[...], "mode":..., "usage":...}
    """
    chosen = agents.route(question, max_specialists=max_specialists)
    route_info = agents.route_explain(question, max_specialists=max_specialists)
    yield {"event": "route", "specialists": chosen,
           "specialist_names": [agents.AGENTS[a]["name"] for a in chosen],
           "scores": route_info["scores"], "dataset_level": route_info["dataset_level"]}

    client = make_client(api_key)

    # ---- Offline fallback: structured multi-specialist report ----
    if client is None:
        yield "\n"
        yield _offline_panel(question, ctx, chosen)
        yield {"event": "final", "trace": [{"id": a, "name": agents.AGENTS[a]["name"]} for a in chosen],
               "mode": "offline", "usage": None}
        return

    # ---- Gather dataset tools (for Data Analyst) ----
    dataset_results = _gather_dataset_results(question, uid) if "data_analyst" in chosen else {}
    dataset_block = _dataset_block(dataset_results)

    # ---- Run specialists in parallel ----
    outputs: list[dict] = []
    with cf.ThreadPoolExecutor(max_workers=max(1, len(chosen))) as ex:
        futs = {ex.submit(_run_specialist, a, question, ctx, dataset_block, client, model): a for a in chosen}
        done_map: dict[str, dict] = {}
        for fut in cf.as_completed(futs):
            o = fut.result()
            done_map[o["id"]] = o
            yield {"event": "specialist_done", "id": o["id"], "name": o["name"],
                   "icon": o["icon"], "ok": bool(o.get("text")) and not o.get("error")}
        # preserve routing order for synthesis
        outputs = [done_map[a] for a in chosen if a in done_map]

    # ---- Synthesize ----
    synth_usage = None
    for piece in _synthesize_stream(question, outputs, ctx, client, model):
        if isinstance(piece, str):
            yield piece
        else:
            synth_usage = piece

    trace = [{"id": o["id"], "name": o["name"], "icon": o["icon"],
              "sources": agents.AGENTS[o["id"]]["sources"],
              "ok": bool(o.get("text")) and not o.get("error"),
              "error": o.get("error"),
              "chars": len(o.get("text", "")),
              "excerpt": (o.get("text", "")[:600])} for o in outputs]
    yield {"event": "final", "trace": trace, "dataset_results": dataset_results,
           "mode": "online", "usage": synth_usage, "model": model}


def _offline_panel(question: str, ctx: dict[str, Any], chosen: list[str]) -> str:
    """No-API fallback: a structured brief from the context + which agents would answer."""
    names = " · ".join(f"{agents.AGENTS[a]['icon']} {agents.AGENTS[a]['name']}" for a in chosen)
    head = (f"**Multi-agent framework (offline mode)**\n\n"
            f"For this question the coordinator would consult: {names}.\n\n"
            f"Set an Anthropic API key to enable the live agent panel. "
            f"Meanwhile, here is the grounded data brief:\n\n---\n\n")
    return head + format_report_card(ctx)


def agents_available() -> dict[str, Any]:
    """UI diagnostic."""
    return {"n_specialists": len(agents.AGENTS),
            "specialists": [{"id": a, **{k: agents.AGENTS[a][k] for k in ("name", "icon", "role")}}
                            for a in agents.SPECIALIST_ORDER],
            "has_key": bool(resolve_api_key())}
