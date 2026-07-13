"""
Chat backend for the protein-page assistant.

Two modes:
- online: Anthropic SDK with streaming. Key resolved from env var
  ANTHROPIC_API_KEY, then st.secrets['anthropic_api_key'].
- offline: no key → render the structured context as a markdown report card.

The public API is `ask(question, ctx, history=None, model=None) -> generator[str, ...]`.
"""
from __future__ import annotations
import os
from typing import Any, Generator, Optional
try:
    from prompts import SYSTEM_PROMPT, SYSTEM_PROMPT_PUBLIC, format_context_prompt
except ImportError:
    from .prompts import SYSTEM_PROMPT, SYSTEM_PROMPT_PUBLIC, format_context_prompt  # type: ignore


DEFAULT_MODEL = "claude-sonnet-4-5"


def resolve_api_key() -> Optional[str]:
    """Try env var, then Streamlit secrets. Returns None if neither is set."""
    key = os.environ.get("ANTHROPIC_API_KEY")
    if key: return key
    try:
        import streamlit as st  # only if running inside streamlit
        key = st.secrets.get("anthropic_api_key") if hasattr(st, "secrets") else None
        if key: return key
    except Exception:
        pass
    return None


def make_client(api_key: Optional[str] = None):
    """Return an anthropic client or None if unavailable."""
    key = api_key or resolve_api_key()
    if not key: return None
    try:
        import anthropic
        return anthropic.Anthropic(api_key=key)
    except ImportError:
        return None


def format_report_card(ctx: dict[str, Any]) -> str:
    """Render the structured context as a markdown report card for offline mode."""
    if "error" in ctx:
        return f"**Cannot look up this protein.** {ctx['error']}"

    ident = ctx["identity"]
    fn = ctx["function"]
    ag = ctx["aging_signal"]
    ml = ctx["ml"]
    orth = ctx["orthologs"]
    doms = ctx["domains"]

    def _yes(v): return v not in (None, "—", "", [])

    lines = []
    lines.append(f"## {ident.get('gene_symbol','?')} ({ident['uniprot_id']}) — {ident.get('recommended_name','?')}")
    lines.append(f"{ident.get('sequence_length','?')} aa · {ident.get('molecular_weight_kda','?')} kDa\n")

    if _yes(fn.get("text")):
        lines.append(f"**Function.** {fn['text']}\n")

    if doms:
        lines.append(f"**Top domains** ({len(doms)} annotated):")
        for d in doms[:5]:
            lines.append(f"- {d['name']} ({d['interpro_id']}, {d['type']}): residues {d['start']}–{d['end']}")
        lines.append("")

    if ag.get("n_significant_aging_peptides"):
        lines.append(f"**Aging signal (LiP-MS).** {ag['n_significant_aging_peptides']} significant peptides "
                     f"(|log2FC|>1, adj-p<0.05, day6+day9) in {len(ag['top_regions'])} contiguous regions:")
        for r in ag["top_regions"]:
            lines.append(f"- Residues {r['start']}–{r['end']} · {r['n_peptides']} peptides · max |log2FC| = {r['max_abs_log2fc']}")
        lines.append("")
    else:
        lines.append("**Aging signal (LiP-MS).** No peptides crossed the significance threshold in aging conditions for this protein.\n")

    if ml.get("top_residues"):
        lines.append(f"**ML top-predicted residues** (model {ml.get('model_version')}, "
                     f"{ml.get('n_high_confidence',0)} with p>0.5):")
        for m in ml["top_residues"][:6]:
            tag = "✓ observed" if m["observed"] else "not observed"
            lines.append(f"- Residue {m['residue']}: p={m['p_destabilized']} · {tag} · observed max |log2FC| = {m['max_abs_log2fc_observed']}")
        lines.append("")

    if orth.get("orthologs"):
        lines.append("**Human orthologs and disease.**")
        for o in orth["orthologs"]:
            lines.append(f"- HGNC:{o['hgnc_id']} ({o['n_methods']} evidence methods)")
        if orth.get("diseases"):
            lines.append(f"\nDisease associations (by orthology): {'; '.join(orth['diseases'][:8])}")
        lines.append("")

    lines.append("---\n*This report was generated offline. Set an Anthropic API key to enable conversational chat about this protein.*")
    return "\n".join(lines)


def ask(question: str, ctx: dict[str, Any],
        history: Optional[list[dict[str, str]]] = None,
        model: str = DEFAULT_MODEL,
        max_tokens: int = 1500,
        api_key: Optional[str] = None,
        audience: str = "researcher") -> Generator[str, None, dict]:
    """
    Stream text chunks for an answer.
    audience: "researcher" (technical, cited [S#]) or "public" (plain language, no bracket codes).
    Yields strings; final yield is a dict summary {'usage': ..., 'model': ..., 'mode': 'online'|'offline'}.
    """
    system_prompt = SYSTEM_PROMPT_PUBLIC if audience == "public" else SYSTEM_PROMPT
    client = make_client(api_key)
    if client is None:
        # Offline: yield the report card as one chunk
        text = format_report_card(ctx)
        yield text
        yield {"usage": None, "model": "offline-report-card", "mode": "offline"}
        return

    user_msg = format_context_prompt(ctx, question)
    messages = list(history or [])
    messages.append({"role": "user", "content": user_msg})

    try:
        with client.messages.stream(
            model=model,
            system=system_prompt,
            max_tokens=max_tokens,
            messages=messages,
        ) as stream:
            for text in stream.text_stream:
                yield text
            final = stream.get_final_message()
            yield {
                "usage": {
                    "input_tokens": final.usage.input_tokens,
                    "output_tokens": final.usage.output_tokens,
                },
                "model": final.model,
                "mode": "online",
                "stop_reason": final.stop_reason,
            }
    except Exception as e:
        yield f"\n\n**Error calling the model:** `{type(e).__name__}: {e}`"
        yield {"usage": None, "model": model, "mode": "error", "error": str(e)}


def api_key_status() -> dict[str, Any]:
    """Diagnostic for the UI to show setup state."""
    key = resolve_api_key()
    source = None
    if os.environ.get("ANTHROPIC_API_KEY"): source = "env var ANTHROPIC_API_KEY"
    elif key:
        try:
            import streamlit as st
            if hasattr(st, "secrets") and st.secrets.get("anthropic_api_key"): source = "streamlit secrets"
        except Exception: pass
    try:
        import anthropic  # noqa
        sdk_ok = True
    except ImportError:
        sdk_ok = False
    return {"has_key": bool(key), "key_source": source, "anthropic_sdk_installed": sdk_ok}
