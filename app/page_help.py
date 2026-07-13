"""
page_help.py — the Help page.

Explains how to use the app (Researcher vs Public), what the colors and heatmap
mean, where the data comes from and how answers are cited, the API-key note for
live chat, and known limitations. Content is drawn from the project READMEs and
the dataset skill so it stays consistent with the science.
"""
from __future__ import annotations
import streamlit as st
import app_core


def _goto(key: str) -> None:
    pages = st.session_state.get("_nav_pages", {})
    target = pages.get(key)
    if target is not None:
        st.switch_page(target)


def render() -> None:
    app_core.top_nav("help")

    st.markdown("## Getting started")
    st.markdown(
        "This app explores a **structural-proteomics map of aging** in the worm *C. elegans*. "
        "A technique called **TMT-LiP-MS** (Sui et al.) measures, for positions along thousands "
        "of proteins, how much each region changes its **shape** between young and old animals "
        "and under related stresses. You reach that data three ways, from the top navigation bar:"
    )
    a, b, c = st.columns(3)
    with a:
        st.markdown("#### 🔬 Researcher")
        st.caption("The full technical viewer + a multi-agent analysis chat with `[S1]–[S8]` "
                   "citations. Search any UniProt ID or gene symbol.")
        if st.button("Open Researcher →", use_container_width=True, key="help_go_res"):
            _goto("researcher")
    with b:
        st.markdown("#### 🌍 Public")
        st.caption("A plain-language, chat-forward tour of aging biology for non-scientists and "
                   "for teaching. Same measurements, jargon-free.")
        if st.button("Explore aging →", use_container_width=True, key="help_go_pub"):
            _goto("public")
    with c:
        st.markdown("#### 🏠 Home")
        st.caption("The overview of what the app is and what it can do.")
        if st.button("Back to Home →", use_container_width=True, key="help_go_home"):
            _goto("home")

    st.markdown("---")

    # ---- Reading the visuals ----
    st.markdown("## Reading the visuals")
    v1, v2 = st.columns(2)
    with v1:
        st.markdown("#### The 3D structure")
        st.markdown(
            "- Each protein is shown as its **AlphaFold** cartoon.\n"
            "- Residues are **painted by log₂ fold-change** in conformational signal:\n"
            "  - 🔴 **red** — shape changes strongly in one direction\n"
            "  - 🔵 **blue** — changes strongly the other way\n"
            "  - ⚪ **grey** — little change, **or not measured**\n"
            "- Spin, zoom, and (in Researcher) focus a region to highlight it in gold."
        )
    with v2:
        st.markdown("#### The 9-condition heatmap")
        st.markdown(
            "Each row is one condition; the x-axis is residue position. Colour is the same "
            "diverging scale as the 3D view. Grouped as:\n"
            "- **Aging** — WT day 6, day 9\n"
            "- **Stress** — heat shock\n"
            "- **PolyQ** — Q35, Q40\n"
            "- **ts-mutants** — myosin / paramyosin at 15 °C and 25 °C\n\n"
            "In Researcher, extra tracks below add an **ML aging-vulnerability** score, "
            "**pathogenic-variant density** (from the human ortholog), and **InterPro domains**."
        )
    st.info("**Grey means *not observed*, not *no change*.** LiP-MS is peptide-level; median "
            "maximum coverage is about 11% of a protein's residues. Absence of colour is absence "
            "of a measurement.")

    st.markdown("---")

    # ---- Where the data comes from + citations ----
    st.markdown("## Where the answers come from")
    st.markdown(
        "Chat answers are **grounded**: the model only reasons over structured context the app "
        "retrieves, and every claim is tagged with a bracketed source. Sources span this study's "
        "measurements and several public databases:"
    )
    st.markdown(
        "| Tag | Source |\n"
        "|---|---|\n"
        "| `[S1]` | UniProt — protein name, function, GO terms, keywords, location |\n"
        "| `[S2]` | InterPro — protein domains and families |\n"
        "| `[S3]` | AlphaFold — predicted 3D structure and pLDDT confidence |\n"
        "| `[S4]` | WormBase SimpleMine — worm→human orthologs |\n"
        "| `[S5]` | **This LiP-MS dataset** — per-residue conformational change in aging/stress |\n"
        "| `[S6]` | ML predictor — per-residue aging-vulnerability probability |\n"
        "| `[S7]` | ClinVar (via UniProt Variation) — disease variants mapped through the ortholog |\n"
        "| `[S8]` | AlphaFold-derived structural context — burial, secondary structure, site distance |\n"
    )
    st.caption("The system prompt explicitly requires the assistant to distinguish **observed** "
               "LiP-MS changes in the worm from **inferred** disease links mapped from the human "
               "ortholog by sequence homology.")

    st.markdown("---")

    # ---- Live chat status (user-facing) ----
    st.markdown("## The conversational chat")
    try:
        import chat_backend
        _ks = chat_backend.api_key_status()
        _has = bool(_ks.get("has_key"))
    except Exception:
        _has = False
    if _has:
        st.markdown("<div class='status-pill live'><span class='dot'></span>"
                    "Live conversational chat is available on the Researcher and Public pages.</div>",
                    unsafe_allow_html=True)
    else:
        st.markdown("<div class='status-pill off'><span class='dot'></span>"
                    "Conversational chat is currently showing an offline grounded summary.</div>",
                    unsafe_allow_html=True)
    st.caption("Either way, the 3D viewer, heatmaps, tracks, and tables are fully available.")

    st.markdown("---")

    # ---- Known limitations ----
    st.markdown("## Known limitations")
    st.markdown(
        "- **Coverage is peptide-level.** Grey residues mean *not observed*; median max coverage "
        "is ~11% per protein.\n"
        "- **Some proteins have no local structure.** 361 measured proteins are TrEMBL entries "
        "newer than the AlphaFold UP000001940 bundle (94.7% of measured proteins are covered).\n"
        "- **Disease links are inferred, not measured.** Variants are mapped from the human "
        "ortholog by pairwise alignment — homology inference, not observations in worms.\n"
        "- **Ortholog inference can mislead for transcription factors.** Always sanity-check TF "
        "orthologs (e.g. SKN-1).\n"
        "- **Disease terms are 'By Orthology'** — they come from the human ortholog's disease "
        "record, not from a worm phenotype."
    )

    st.markdown("---")
    st.markdown("## Good first questions")
    st.markdown(
        "- *What is this protein and what does it do?*\n"
        "- *Which regions change with age, and where do they sit in the domain structure?*\n"
        "- *Is this change aging-specific, or also seen under heat shock and polyQ?*\n"
        "- *Is the aging region buried in the fold, and could that affect stability?*\n"
        "- *Are there human disease variants near the aging region?*"
    )
    st.caption("Dataset: 6,823 proteins × 9 conditions from TMT-LiP-MS in *C. elegans* (Sui et al.). "
               "Structures from AlphaFold UP000001940.")
