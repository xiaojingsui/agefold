"""
public_view.py — renders the public "Explore aging" track.

Reuses the same data loaders and structure-rendering approach as the researcher
view, but re-narrates everything in plain language. Rendered by the Public page
(page_public.py) in the top-nav shell.

The caller passes in the shared ctx (built by app_core.build_ctx) so we don't
duplicate caching logic. See render_public() signature.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
import streamlit as st
import py3Dmol
import streamlit.components.v1 as components

import public_content


def render_public(ctx: dict) -> None:
    """Render the entire public track.

    ctx carries shared resources from app_core.build_ctx():
      prot_feat, gene2uid, per_res, af_index (dict), read_pdb (fn),
      diverging_color (fn), load_variants (fn), load_orthologs (fn),
      CONDITIONS, COND_LABEL
    """
    # Router: landing page vs a single featured/selected protein.
    sel = st.session_state.get("public_uid")
    if sel:
        _render_protein(ctx, sel)
    else:
        _render_landing(ctx)


def _render_landing(ctx: dict) -> None:
    from pathlib import Path
    L = public_content.LANDING

    # ---- Chat-forward entry: lead with a conversation, not a wall of text ----
    st.markdown(
        "<div class='hero-head'><h2 class='hero-title'>Chat your way into the biology of aging</h2>"
        "<p class='hero-sub'>Type the name of a protein or gene and start a plain-language "
        "conversation — what it does, how its shape changes as a worm grows old, and what that "
        "might mean for human health. No science background needed.</p></div>",
        unsafe_allow_html=True)

    prot_feat = ctx["prot_feat"]
    gene2uid = ctx["gene2uid"]
    _uidset = set(prot_feat["uniprot_id"].values)

    with st.form("public_start_chat", clear_on_submit=False):
        c_in, c_btn = st.columns([4, 1])
        with c_in:
            _q = st.text_input(
                "Protein or gene name",
                label_visibility="collapsed",
                placeholder="Try a gene name like  sod-1  ·  hsp-1  ·  daf-16  ·  vit-6",
            )
        with c_btn:
            _go = st.form_submit_button("Start chat →", use_container_width=True, type="primary")
    if _go and _q:
        q = _q.strip()
        hit = None
        if q.upper() in _uidset:
            hit = q.upper()
        elif q.lower() in gene2uid:
            hit = gene2uid[q.lower()]
        if hit:
            st.session_state["public_uid"] = hit
            st.rerun()
        else:
            st.warning(f"I couldn't find “{q}”. Try a gene name like **sod-1**, **hsp-1**, or **vit-6**, "
                       "or pick one of the featured stories below.")

    # Quick-start chips → jump straight into a conversation
    st.caption("Or jump straight in:")
    _chip_cols = st.columns(4)
    for _i, (_uid, _gene, _blurb) in enumerate(public_content.QUICKSTART if hasattr(public_content, "QUICKSTART")
                                               else []):
        if _chip_cols[_i % 4].button(f"💬 {_gene}", key=f"qs_{_uid}", use_container_width=True,
                                     help=_blurb):
            st.session_state["public_uid"] = _uid
            st.rerun()

    st.markdown("---")

    st.markdown(f"### 🌍 {L['headline']}")
    st.markdown(L["intro"])
    st.caption(f"📊 {L['dataset_line']}")

    st.markdown("---")

    # How to read
    st.markdown("### How to read what you'll see")
    diagram = Path(__file__).parent / "data" / "public_howto_diagram.png"
    if not diagram.exists():
        diagram = Path(__file__).parent.parent / "data" / "public_howto_diagram.png"
    if diagram.exists():
        # Center at 1/3 width so the diagram's height is ~1/3 of full-width.
        _dl, _dc, _dr = st.columns([1, 1, 1])
        with _dc:
            st.image(str(diagram), use_container_width=True)
    st.markdown(
        "Scientists measured, for **every position** along each protein, how much its shape "
        "shifts between young and old worms. Bright colors mark the parts that change the most."
    )

    # Why worm
    with st.container():
        st.info(L["why_worm"])

    st.markdown("---")

    # Featured gallery
    st.markdown("### ✨ Featured stories — click any card to explore")
    feats = public_content.FEATURED
    # 2-column grid of cards
    for i in range(0, len(feats), 2):
        cols = st.columns(2)
        for j, col in enumerate(cols):
            if i + j >= len(feats):
                continue
            f = feats[i + j]
            with col:
                with st.container(border=True):
                    st.markdown(f"#### {f['icon']} {f['title']}")
                    st.markdown(f"*{f['hook']}*")
                    if "human_health" in f:
                        st.caption("🧬 Linked to a human disease")
                    if st.button(f"Explore {f['gene']} →", key=f"feat_{f['uid']}",
                                 use_container_width=True):
                        st.session_state["public_uid"] = f["uid"]
                        st.rerun()

    st.markdown("---")

    # Glossary
    with st.expander("📖 Aging science basics — plain-language glossary"):
        for term, definition in public_content.GLOSSARY:
            st.markdown(f"**{term}** — {definition}")

    st.caption("This is the public 'Explore aging' view. Researchers can open the full "
               "technical view from the **Researcher** tab in the top navigation bar.")


def _featured_entry(uid: str) -> dict | None:
    for f in public_content.FEATURED:
        if f["uid"] == uid:
            return f
    return None


def _plain_regions(per_res: pd.DataFrame, uid: str, seq_len: int, top_k: int = 4) -> list[tuple[int, int]]:
    """Contiguous aging-changed spans (day6/day9, |log2FC|>1, adj-p<0.05),
    merged within 20 aa, top-k by peak magnitude — for plain-language narration."""
    a = per_res[(per_res["uniprot_id"] == uid)
                & (per_res["condition"].isin(["day6", "day9"]))
                & (per_res["log2fc_max_abs"].abs() > 1.0)
                & (per_res["adj_pval_at_max_abs"] < 0.05)
                & (per_res["adj_pval_at_max_abs"] >= 0)]
    if not len(a):
        return []
    resis = sorted(a["residue"].unique())
    spans = []
    start = prev = resis[0]
    for r in resis[1:]:
        if r - prev <= 20:
            prev = r
        else:
            spans.append((start, prev)); start = prev = r
    spans.append((start, prev))
    # rank by peak |log2FC| within span
    def peak(s):
        seg = a[(a["residue"] >= s[0]) & (a["residue"] <= s[1])]
        return seg["log2fc_max_abs"].abs().max()
    spans.sort(key=peak, reverse=True)
    return spans[:top_k]


def _render_protein(ctx: dict, uid: str) -> None:
    prot_feat = ctx["prot_feat"]
    per_res = ctx["per_res"]
    af_index = ctx["af_index"]

    if st.button("← Back to featured proteins"):
        st.session_state["public_uid"] = None
        st.session_state.pop("public_pick", None)
        st.rerun()

    row = prot_feat[prot_feat["uniprot_id"] == uid]
    if not len(row):
        st.error("That protein isn't in the dataset.")
        return
    prow = row.iloc[0]
    gene = prow["gene_symbol"] if pd.notna(prow.get("gene_symbol")) else uid
    seq_len = int(prow["seq_length"]) if pd.notna(prow.get("seq_length")) else None

    feat = _featured_entry(uid)

    # Header
    title = feat["title"] if feat else f"{gene}"
    icon = feat["icon"] if feat else "🔬"
    st.markdown(f"# {icon} {title}")
    if feat:
        st.markdown(f"### *{feat['hook']}*")

    # Two columns: structure | story
    col_struct, col_story = st.columns([1.3, 1])

    with col_struct:
        st.markdown("#### The 3D shape — spin it around")
        pdb_text = ctx["get_pdb_text"](uid)
        regions = _plain_regions(per_res, uid, seq_len) if seq_len else []
        if pdb_text:
            # color by day6 aging signal
            day6 = per_res[(per_res["uniprot_id"] == uid) & (per_res["condition"] == "day6")]
            vals = day6["log2fc_max_abs"].dropna().values
            vmax = max(float(np.nanpercentile(np.abs(vals), 98)) if len(vals) else 1.0, 0.5)
            res2color = {int(r): ctx["diverging_color"](v, vmax)
                         for r, v in zip(day6["residue"], day6["log2fc_max_abs"]) if np.isfinite(v)}
            view = py3Dmol.view(width=560, height=460)
            view.addModel(pdb_text, "pdb")
            view.setStyle({"cartoon": {"color": "#dcdcdc"}})
            for r, c in res2color.items():
                view.setStyle({"resi": str(r)}, {"cartoon": {"color": c}})
            view.zoomTo(); view.setBackgroundColor("white")
            components.html(view._make_html(), height=480, scrolling=False)
            st.caption("🔴 Red = this part gets **tighter / more protected** with age   ·   "
                       "🔵 Blue = this part gets **looser / more open** with age   ·   "
                       "⚪ Grey = barely changes, or not measured")
        else:
            st.warning("No 3D structure is available for this protein.")

    with col_story:
        # What it does — plain language, from featured story or UniProt function
        st.markdown("#### What this protein does")
        if feat:
            st.markdown(feat["story"])
        else:
            fn = _uniprot_function_plain(ctx, uid)
            st.markdown(fn or "_A description for this protein isn't available yet._")

        # Where it changes
        st.markdown("#### Where it changes as the worm ages")
        if regions:
            span_txt = ", ".join(f"**{s}–{e}**" for s, e in regions)
            st.markdown(
                f"The biggest age-related shape changes happen in these stretches of the protein: "
                f"{span_txt}. On the 3D structure, look for the brightly colored patches."
            )
        else:
            st.markdown("_This protein shows little measurable shape change with age in our data._")

    # Full-width sections below
    if feat and "human_health" in feat:
        st.markdown("---")
        st.markdown("#### 🧬 The human health connection")
        st.markdown(feat["human_health"])
        st.caption("Based on matching the worm protein to its human counterpart and known disease "
                   "mutations. This is an inference from similarity, not a measurement in humans.")

    if feat:
        st.markdown("---")
        st.markdown("#### 💡 Why this matters for aging")
        st.success(feat["why_matters"])

    # Public chat
    st.markdown("---")
    _render_public_chat(ctx, uid, gene)

    # Researcher escape hatch
    st.markdown("---")
    with st.expander("🔬 For researchers — see the technical view"):
        st.markdown(
            f"To open **{gene} ({uid})** in the full technical viewer — with the 9-condition "
            f"heatmap, the machine-learning aging-vulnerability track, mapped disease variants, "
            f"and protein-domain annotations — open the **🔬 Researcher** tab in the top "
            f"navigation bar, then search for `{uid}`."
        )
        st.markdown("**Technical summary (this protein):**")
        peak_regions = _plain_regions(per_res, uid, seq_len) if seq_len else []
        if peak_regions:
            for s, e in peak_regions:
                seg = per_res[(per_res["uniprot_id"] == uid)
                              & (per_res["condition"].isin(["day6", "day9"]))
                              & (per_res["residue"] >= s) & (per_res["residue"] <= e)]
                peak = seg["log2fc_max_abs"].abs().max()
                st.markdown(f"- residues **{s}–{e}**: peak |log2FC| = {peak:.2f} (day6/day9, adj-p < 0.05)")
        else:
            st.markdown("- No significant aging regions (|log2FC| > 1, adj-p < 0.05).")

    # Glossary always available
    with st.expander("📖 Aging science basics — plain-language glossary"):
        for term, definition in public_content.GLOSSARY:
            st.markdown(f"**{term}** — {definition}")


def _render_public_chat(ctx: dict, uid: str, gene: str) -> None:
    st.markdown("#### 💬 Ask a question about this protein")
    st.caption("Ask anything — the answer is written for curious non-scientists and is based on "
               "the same data researchers see.")
    try:
        import retrieval, chat_backend
    except ImportError as e:
        st.warning(f"Chat isn't available ({e}).")
        return

    pctx = retrieval.get_protein_context(uid)
    key = f"public_chat::{uid}"
    if key not in st.session_state:
        st.session_state[key] = []

    # Suggestion chips
    suggestions = [
        f"What does {gene} do, in simple terms?",
        "Why does this protein matter for aging?",
        "Where does it change as the worm gets older?",
    ]
    if any("human_health" in f for f in public_content.FEATURED if f["uid"] == uid):
        suggestions.append("Is this linked to any human disease?")

    cols = st.columns(len(suggestions))
    clicked = None
    for c, s in zip(cols, suggestions):
        if c.button(s, key=f"pub_sugg_{uid}_{hash(s)%9999}", use_container_width=True):
            clicked = s

    # Replay history
    for msg in st.session_state[key]:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    user_input = st.chat_input("Type your question…", key=f"pub_chatin_{uid}") or clicked
    if user_input:
        st.session_state[key].append({"role": "user", "content": user_input})
        with st.chat_message("user"):
            st.markdown(user_input)
        with st.chat_message("assistant"):
            placeholder = st.empty()
            acc = ""
            key_status = chat_backend.api_key_status()
            if key_status["has_key"]:
                try:
                    gen = chat_backend.ask(
                        user_input, pctx,
                        history=[{"role": m["role"], "content": m["content"]}
                                 for m in st.session_state[key][:-1]],
                        max_tokens=1200, audience="public",
                    )
                    for chunk in gen:
                        if isinstance(chunk, dict):
                            break
                        acc += chunk
                        placeholder.markdown(acc + "▌")
                    placeholder.markdown(acc)
                except Exception as e:
                    acc = f"Sorry — something went wrong: {e}"
                    placeholder.markdown(acc)
            else:
                acc = ("💡 The conversational answers need an API key, which isn't set up here. "
                       "But you can still explore this protein's story and 3D shape above!")
                placeholder.markdown(acc)
        st.session_state[key].append({"role": "assistant", "content": acc})

    # Sources under the hood
    with st.expander("Where does this information come from?"):
        st.markdown(
            "Answers draw on public protein databases (UniProt, InterPro, AlphaFold), the "
            "WormBase worm genetics resource, this study's LiP-MS aging measurements, and — where "
            "a human disease link exists — ClinVar variant records matched through the human version "
            "of the gene. Nothing is made up; if the data doesn't cover something, the assistant says so."
        )


def _uniprot_function_plain(ctx: dict, uid: str) -> str:
    """Best-effort plain function text from UniProt cache (non-featured proteins)."""
    try:
        u = ctx["cached_uniprot"](uid)
        txt = u.get("function_text") or u.get("function") or ""
        # Trim citation clutter like (PubMed:12345)
        import re
        txt = re.sub(r"\(PubMed:[^)]*\)", "", txt)
        txt = re.sub(r"\(By similarity\)", "", txt)
        txt = re.sub(r"\s+", " ", txt).strip()
        return txt[:400] + ("…" if len(txt) > 400 else "")
    except Exception:
        return ""
