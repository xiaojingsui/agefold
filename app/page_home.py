"""
page_home.py — the Home showcase / landing page.

Leads with the product promise ("A conversational AI to decode the biology of
aging"), then dataset stat pills and three capability cards, each with a call to
action that jumps to the relevant page. Cross-page jumps use the st.Page objects
the shell stashes in st.session_state["_nav_pages"].
"""
from __future__ import annotations
import base64
from pathlib import Path
import streamlit as st
import app_core


@st.cache_data
def _crane_data_uri() -> str:
    """Base64-encode the origami-crane folding sequence for inline hero use."""
    p = Path(__file__).parent / "assets" / "crane_folding_steps.png"
    if not p.exists():
        return ""
    b64 = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def _goto(key: str) -> None:
    pages = st.session_state.get("_nav_pages", {})
    target = pages.get(key)
    if target is not None:
        st.switch_page(target)


def render() -> None:
    app_core.top_nav("home")
    ctx = app_core.build_ctx()
    per_res = ctx["per_res"]

    # ---- Hero: conversational-AI-first framing ----
    _crane = _crane_data_uri()
    _crane_img = (
        f'<img class="hero-crane" src="{_crane}" '
        'alt="Folding an origami crane, step by step">' if _crane else ""
    )
    st.markdown(
        f"""
<div class="land-hero">
  <p class="eyebrow">C. elegans · structural proteomics · aging</p>
  <h1>A <span class="accent">conversational AI</span> to decode the biology of aging</h1>
  <p class="lede">Ask questions in plain language and get grounded, cited answers about how
  thousands of proteins change shape as an animal grows old — measured by TMT-LiP-MS across
  nine conditions and painted onto 3D AlphaFold structures. Built for experts and the
  curious public alike.</p>
  {_crane_img}
</div>
""",
        unsafe_allow_html=True,
    )

    # ---- Primary CTAs ----
    c1, c2, c3 = st.columns([1, 1, 1])
    with c1:
        if st.button("🔬  Open the Researcher workspace", use_container_width=True):
            _goto("researcher")
    with c2:
        if st.button("🌍  Explore aging (for everyone)", use_container_width=True):
            _goto("public")
    with c3:
        if st.button("❓  How it works", use_container_width=True):
            _goto("help")

    # ---- Dataset stat pills ----
    st.markdown("<div style='height:10px'></div>", unsafe_allow_html=True)
    try:
        n_measured = int(per_res["uniprot_id"].nunique())
    except Exception:
        n_measured = 6823
    try:
        import agents as _agents
        n_agents = len(_agents.AGENTS)
    except Exception:
        n_agents = 4
    ctx["stat_pills"]([
        (f"{n_measured:,}", "proteins measured"),
        ("9", "conditions"),
        ("530,947", "peptide measurements"),
        ("94.7%", "with AlphaFold structure"),
        ("58%", "with human ortholog"),
        (f"{n_agents}+", "specialist AI agents"),
    ])

    st.markdown("---")

    # ---- What it can do: three capability cards ----
    st.markdown("### What you can do here")
    st.caption("One dataset, three ways in — all reading the same measurements.")

    fc1, fc2, fc3 = st.columns(3)
    with fc1:
        st.markdown(
            """
<div class="feat-card">
  <div class="fc-icon">🔬</div>
  <div class="fc-title">Talk to a specialist panel</div>
  <div class="fc-body">In <b>Researcher</b> mode, a coordinator routes your question to
  structural-biology, aging &amp; stress, disease-genetics and data-analyst agents that
  reason over the grounded data and return one answer with <code>[S1]–[S8]</code> citations.</div>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("Open Researcher →", key="card_researcher", use_container_width=True):
            _goto("researcher")
    with fc2:
        st.markdown(
            """
<div class="feat-card">
  <div class="fc-icon">🌍</div>
  <div class="fc-title">Explore aging in plain language</div>
  <div class="fc-body">In <b>Public</b> mode, anyone can chat their way into the science —
  what a protein does, where its shape changes with age, and the human-health connection.
  Written for curious non-scientists and great for teaching.</div>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("Explore aging →", key="card_public", use_container_width=True):
            _goto("public")
    with fc3:
        st.markdown(
            """
<div class="feat-card">
  <div class="fc-icon">🧬</div>
  <div class="fc-title">See it on the 3D structure</div>
  <div class="fc-body">Every answer is anchored in an interactive <b>data viewer</b>: per-residue
  conformational change painted on AlphaFold cartoons, a 9-condition heatmap, InterPro
  domains, an ML aging-vulnerability track, and mapped disease variants.</div>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("Open the viewer →", key="card_viewer", use_container_width=True):
            _goto("researcher")

    st.markdown("---")

    # ---- How it works, three steps ----
    st.markdown("### How it works")
    s1, s2, s3 = st.columns(3)
    with s1:
        st.markdown("**1 · Ask**")
        st.caption("Type a protein and a question — in expert or plain language.")
    with s2:
        st.markdown("**2 · Grounded retrieval**")
        st.caption("The app pulls this study's LiP-MS measurements plus live UniProt, InterPro, "
                   "AlphaFold, WormBase and ClinVar context.")
    with s3:
        st.markdown("**3 · Cited answer + 3D**")
        st.caption("You get a synthesized, source-cited answer alongside the region lit up on "
                   "the 3D structure.")

    st.markdown(
        "<div style='height:6px'></div>", unsafe_allow_html=True)
    st.caption("Data: TMT-LiP-MS conformational proteomics in *C. elegans* (Sui et al.). "
               "Structures from AlphaFold UP000001940. See the **Help** tab for provenance, "
               "the meaning of the colors, and known limitations.")
