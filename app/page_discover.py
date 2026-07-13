"""
page_discover.py — the "Discover" page.

A full-width home for the composite candidate ranking: proteins scored by their
aging conformational signal combined with human disease-variant overlap. Filters,
a ranked table, and a jump box that loads any hit straight into the Researcher
viewer. Pulled out of the Researcher page so discovery is its own destination.
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
    app_core.top_nav("discover")
    ctx = app_core.build_ctx()
    load_discovery = ctx["load_discovery"]

    st.markdown("### 🔍 Discovery — top-ranked candidates")
    st.caption("Composite score combines aging conformational signal + human disease-variant overlap. "
                "See app/DISCOVERY_README.md for the scoring formula.")

    disc = load_discovery()
    if not len(disc):
        st.warning("No discovery scores found. Run `python app/discovery.py` to generate them.")
    else:
        with st.container():
            cf1, cf2, cf3 = st.columns([1.2, 1.2, 1])
            min_ident = cf1.slider("Min alignment identity (worm↔human)",
                                    min_value=0.0, max_value=1.0, value=0.3, step=0.05,
                                    help="Filters out low-confidence variant mappings. 0.3 is the recommended floor.")
            min_composite = cf2.slider("Min composite score",
                                        min_value=0.0, max_value=float(disc["composite"].max()),
                                        value=0.0, step=1.0)
            top_n = cf3.number_input("Show top N", min_value=10, max_value=500, value=50, step=10)

        cf4, cf5, cf6 = st.columns([1.5, 1.5, 1])
        require_variant = cf4.checkbox("Require ≥1 pathogenic variant mapped", value=False,
                                        help="Restricts to proteins whose human ortholog has ClinVar pathogenic entries.")
        require_overlap = cf5.checkbox("Require aging-variant overlap ≥1", value=False,
                                        help="Restricts to proteins where pathogenic variants fall inside aging regions.")

        d = disc.copy()
        d = d[(d["alignment_identity"].fillna(0) >= min_ident) & (d["composite"] >= min_composite)]
        if require_variant:
            d = d[d["variant_load"] > 0]
        if require_overlap:
            d = d[d["aging_variant_overlap"] > 0]
        st.caption(f"**{len(d):,} proteins pass filters** (of {len(disc):,} total). Showing top {min(top_n, len(d))}.")

        show = d.head(top_n).copy()
        show["composite"] = show["composite"].round(2)
        show["hotspot_max"] = show["aging_hotspot_max"].round(2)
        show["load_per_aa"] = show["aging_load_per_aa"].round(3)
        show["alignment_id"] = show["alignment_identity"].round(2)
        show["disease"] = show["top_disease"].astype(str).where(show["top_disease"].notna(), "")
        show_cols = ["uniprot_id", "gene_symbol", "seq_length", "composite",
                      "hotspot_max", "load_per_aa", "variant_load", "aging_variant_overlap",
                      "alignment_id", "disease"]
        show["gene"] = show["gene_symbol"].astype(str).where(show["gene_symbol"].notna(), "")
        display = show[show_cols].rename(columns={
            "uniprot_id": "UniProt", "gene_symbol": "gene", "seq_length": "aa",
            "hotspot_max": "hotspot", "load_per_aa": "load/aa",
            "variant_load": "path.var.", "aging_variant_overlap": "aging∩var",
            "alignment_id": "aln id",
        })
        st.dataframe(display, hide_index=True, use_container_width=True, height=520)

        st.markdown("**Drill down:** copy a UniProt ID from the first column into the sidebar to load that protein's page. "
                    "Or paste one below to jump directly:")
        jump_col, _ = st.columns([2, 3])
        jump = jump_col.text_input("Jump to UniProt ID from the table", value="",
                                    placeholder="e.g. P34697", key="_jump_uid")
        if jump:
            j = jump.strip().upper()
            if j in set(disc["uniprot_id"].values):
                st.session_state["protein_query"] = j
                _goto("researcher")
                st.rerun()
            else:
                st.error(f"'{jump}' is not in the scored table.")
