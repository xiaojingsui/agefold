"""
researcher_view.py — the full technical viewer, refactored out of the old
streamlit_app.py monolith into a page module.

render(ctx) draws a single vertical research flow:
  • sidebar: protein search, condition selector, painting controls
  • Overview header, then a two-column row: Protein card (name, function,
    location) on the left, UniProt information (IDs, gene, length, MW, GO,
    keywords) on the right
  • the 3D AlphaFold structure (with region-focus controls + colorbar) and the
    multi-track per-residue heatmap (9 conditions + ML / variant / InterPro-
    domain tracks), full width below the card
  • the structural-context readout below the visuals
  • ML top residues, disease variants, and the human-ortholog panel
  • Conversation: one unified multi-agent chat that routes to specialist agents
    and returns cited answers
  (Discover, the composite candidate ranking, is its own top-nav page:
  page_discover.py)

Every data resource comes from ctx (built in app_core.build_ctx); no module-level
data globals. Session-state keys are preserved so the focus selection, the
conversation history, and the Discover→viewer jump keep working.
"""
from __future__ import annotations
import os, re
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
import py3Dmol
import streamlit.components.v1 as components


def render(ctx: dict) -> None:
    # ---- unpack shared resources ----
    prot_feat = ctx["prot_feat"]
    gene2uid = ctx["gene2uid"]
    per_res = ctx["per_res"]
    interpro = ctx["interpro"]
    af_index = ctx["af_index"]
    read_pdb = ctx["read_pdb"]
    get_pdb_text = ctx["get_pdb_text"]
    diverging_color = ctx["diverging_color"]
    load_variants = ctx["load_variants"]
    load_orthologs = ctx["load_orthologs"]
    load_predictions = ctx["load_predictions"]
    cached_uniprot = ctx["cached_uniprot"]
    cached_hgnc = ctx["cached_hgnc"]
    CONDITIONS = ctx["CONDITIONS"]
    COND_LABEL = ctx["COND_LABEL"]
    COND_GROUP = ctx["COND_GROUP"]
    LIPMS_CMAP = ctx["LIPMS_CMAP"]
    _stat_pills = ctx["stat_pills"]

    # =================== Sidebar (protein & painting) ===================
    st.sidebar.markdown(
        "<div class='sb-brand'><div class='sb-mark'>🧬</div>"
        "<div><div class='sb-name'>LiP-MS Aging Viewer</div>"
        "<div class='sb-tag'>C. elegans · structural proteomics</div></div></div>",
        unsafe_allow_html=True)
    st.sidebar.markdown("<div class='sb-label'>Protein &amp; painting</div>", unsafe_allow_html=True)
    st.sidebar.caption("Search a protein, choose a condition, and tune how the structure is painted.")

    _default_uid = st.session_state.get("protein_query", "P18948") or "P18948"
    query = st.sidebar.text_input(
        "Protein (UniProt ID or gene symbol)",
        value=_default_uid,
        help="Try P18948 (vit-6), P09446 (hsp-1), P41988 (skn-1), P52810 (act-1). "
             "Or use the Discover page's jump box to load a top hit.",
    )
    st.session_state["protein_query"] = query

    uid = None
    if query:
        q = query.strip()
        q_upper = q.upper()
        if q_upper in set(prot_feat["uniprot_id"].values):
            uid = q_upper
        elif q.lower() in gene2uid:
            uid = gene2uid[q.lower()]

    if uid is None:
        st.sidebar.error(f"'{query}' not found. Try a UniProt accession or gene symbol.")
        st.stop()

    cond = st.sidebar.selectbox(
        "Condition",
        options=[c for c, _, _ in CONDITIONS],
        format_func=lambda c: f"{COND_LABEL[c]}  ({COND_GROUP[c]})",
        index=0,
    )

    st.sidebar.markdown("<div class='sb-label'>Painting</div>", unsafe_allow_html=True)
    col_metric = st.sidebar.radio("Residue value", ["log2FC max |·|", "log2FC mean"], horizontal=True, index=0)
    col_metric_col = "log2fc_max_abs" if col_metric.startswith("log2FC max") else "log2fc_mean"

    auto_scale = st.sidebar.checkbox("Auto scale (98th percentile)", value=True)
    vmax_input = None if auto_scale else st.sidebar.slider("|log2FC| max", min_value=0.5, max_value=6.0, value=3.0, step=0.1)
    sig_only = st.sidebar.checkbox("Grey out non-significant (adj-p ≥ 0.05)", value=False)

    # ---------------- Data slice ----------------
    prow = prot_feat[prot_feat["uniprot_id"] == uid].iloc[0]
    gene = prow["gene_symbol"] if pd.notna(prow["gene_symbol"]) else uid
    seq_len = int(prow["seq_length"]) if pd.notna(prow["seq_length"]) else None

    sub = per_res[(per_res["uniprot_id"] == uid) & (per_res["condition"] == cond)].copy()
    if sig_only and len(sub):
        sub.loc[sub["adj_pval_at_max_abs"] >= 0.05, col_metric_col] = np.nan

    # ---------------- Header ----------------
    _flags = []
    if pd.notna(prow.get("longevity")):
        _flags.append(f"🧓 {prow['longevity']}")
    if pd.notna(prow.get("hsp90_client")):
        _flags.append("🔗 HSP90 client")
    if pd.notna(prow.get("tric_client")):
        _flags.append("🔗 TRiC client")
    if _flags:
        st.markdown("<div class='flag-row'>" + "".join(f"<span class='chip'>{f}</span>" for f in _flags) + "</div>",
                    unsafe_allow_html=True)


    # =============== local helpers (close over ctx data) ===============
    _UIDSET = set(prot_feat["uniprot_id"].values)

    def _aging_regions_for(puid, merge_gap=20, top_n=6):
        a = per_res[(per_res["uniprot_id"] == puid)
                    & (per_res["condition"].isin(["day6", "day9"]))
                    & (per_res["log2fc_max_abs"].abs() > 1.0)
                    & (per_res["adj_pval_at_max_abs"] < 0.05)
                    & (per_res["adj_pval_at_max_abs"] >= 0)]
        resis = sorted(a["residue"].astype(int).unique())
        if not resis:
            return []
        regs = []
        s = prev = resis[0]
        for r in resis[1:]:
            if r - prev <= merge_gap:
                prev = r
            else:
                regs.append((s, prev)); s = prev = r
        regs.append((s, prev))
        def _peak(rg):
            seg = a[(a["residue"] >= rg[0]) & (a["residue"] <= rg[1])]
            return float(seg["log2fc_max_abs"].abs().max())
        return sorted(regs, key=_peak, reverse=True)[:top_n]

    def _parse_selection(text):
        if not text:
            return None
        t = text.strip().replace("–", "-").replace(" ", "")
        try:
            if "-" in t:
                a, b = t.split("-", 1)
                a, b = int(a), int(b)
                return (min(a, b), max(a, b))
            v = int(t)
            return (v, v)
        except ValueError:
            return None

    def _resolve_uid_from_text(text, fallback=None):
        if not text:
            return fallback
        _acc = r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z0-9]{3}[0-9]){1,2})\b"
        for tok in re.findall(_acc, text.upper()):
            if tok in _UIDSET:
                return tok
        for w in re.findall(r"[A-Za-z0-9\-\.]+", text.lower()):
            if w in gene2uid:
                return gene2uid[w]
        return fallback

    def _render_3d_panel(puid, focus=None, height=360):
        pdb_text = get_pdb_text(puid)
        if not pdb_text:
            st.info(f"No AlphaFold structure available for {puid}.")
            return
        d = per_res[(per_res["uniprot_id"] == puid) & (per_res["condition"] == "day6")]
        vals = d["log2fc_max_abs"].dropna().values
        vmax = max(float(np.nanpercentile(np.abs(vals), 98)) if len(vals) else 1.0, 0.5)
        res2color = {int(r): diverging_color(v, vmax)
                     for r, v in zip(d["residue"], d["log2fc_max_abs"]) if np.isfinite(v)}
        view = py3Dmol.view(width=430, height=height)
        view.addModel(pdb_text, "pdb")
        view.setStyle({"cartoon": {"color": "#dcdcdc"}})
        for r, c in res2color.items():
            view.setStyle({"resi": str(r)}, {"cartoon": {"color": c}})
        if focus:
            rr = f"{int(focus[0])}-{int(focus[1])}"
            view.addStyle({"resi": rr}, {"cartoon": {"thickness": 1.0},
                                         "stick": {"radius": 0.18, "color": "#e0a300"}})
            view.zoomTo({"resi": rr})
        else:
            view.zoomTo()
        view.setBackgroundColor("white")
        components.html(view._make_html(), height=height + 12, scrolling=False)

    def _render_heatmap_panel(puid, focus=None):
        pl = prot_feat[prot_feat["uniprot_id"] == puid]
        if not len(pl):
            return
        slen = int(pl.iloc[0]["seq_length"]) if pd.notna(pl.iloc[0]["seq_length"]) else None
        if not slen:
            return
        cond_order = [c for c, _, _ in CONDITIONS]
        mat = np.full((len(cond_order), slen), np.nan)
        pa = per_res[per_res["uniprot_id"] == puid]
        for i, c in enumerate(cond_order):
            s = pa[pa["condition"] == c]
            if len(s):
                m = (s["residue"].values >= 1) & (s["residue"].values <= slen)
                mat[i, s["residue"].values[m] - 1] = s["log2fc_max_abs"].values[m]
        finite = mat[np.isfinite(mat)]
        vmax = max(float(np.nanpercentile(np.abs(finite), 98)) if finite.size else 1.0, 0.5)
        fig, ax = plt.subplots(figsize=(7, 2.9))
        norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
        ax.imshow(mat, aspect="auto", cmap=LIPMS_CMAP, norm=norm, interpolation="nearest",
                  extent=[0.5, slen + 0.5, len(cond_order) - 0.5, -0.5])
        ax.set_yticks(range(len(cond_order)))
        ax.set_yticklabels([COND_LABEL[c] for c in cond_order], fontsize=7)
        ax.set_xlabel("Residue position", fontsize=8)
        if focus:
            ax.axvspan(focus[0] - 0.5, focus[1] + 0.5, facecolor="none", edgecolor="#e0a300", linewidth=2)
        for spn in ax.spines.values():
            spn.set_edgecolor("#bbbbbb")
        fig.tight_layout()
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)

    # ================= shared visualization state + two-column layout =================
    _aging_regions = _aging_regions_for(uid)
    _sel_key = f"focus_resi::{uid}"
    focus_sel = st.session_state.get(_sel_key)

    # ===================== Overview: card + UniProt (top row) =====================
    st.markdown("## Overview")
    uni = cached_uniprot(uid)
    _card_col, _uni_col = st.columns([1, 1], gap="large")
    with _card_col:
        st.markdown("### Protein card")
        if uni.get("recommended_name"):
            st.markdown(f"**{uni['recommended_name']}**")
        if uni.get("function"):
            fn = uni["function"]
            st.markdown(f"<div style='font-size:0.9em'>{fn[:400]}{'…' if len(fn)>400 else ''}</div>",
                        unsafe_allow_html=True)
        if uni.get("subcellular_location"):
            st.caption("📍 " + ", ".join(uni["subcellular_location"][:4]))
        if pd.notna(prow.get("tissue_specificity")) and prow["tissue_specificity"] != "not tissue-specific":
            st.markdown(f"<div class='flag-row'><span class='chip'>🩸 {prow['tissue_specificity']}</span></div>",
                        unsafe_allow_html=True)
    with _uni_col:
        st.markdown("### UniProt information")
        st.markdown(f"**UniProt**: [`{uid}`](https://www.uniprot.org/uniprotkb/{uid})")
        st.markdown(f"**Gene symbol**: *{gene}*")
        if seq_len:
            st.markdown(f"**Sequence length**: {seq_len:,} aa")
        if pd.notna(prow.get("mol_weight")):
            st.markdown(f"**MW**: {prow['mol_weight']/1000:.1f} kDa")
        if uni.get("go_bp") or uni.get("go_mf") or uni.get("go_cc"):
            with st.expander(f"GO terms ({len(uni['go_bp'])+len(uni['go_mf'])+len(uni['go_cc'])})", expanded=False):
                for ns_key, ns_name in [("go_bp", "Biological Process"),
                                         ("go_mf", "Molecular Function"),
                                         ("go_cc", "Cellular Component")]:
                    terms = uni.get(ns_key, [])
                    if terms:
                        st.markdown(f"**{ns_name}**")
                        for t in terms[:8]:
                            st.markdown(f"- [{t['id']}](https://www.ebi.ac.uk/QuickGO/term/{t['id']}) — {t['name']}")
        if uni.get("keywords"):
            with st.expander(f"UniProt keywords ({len(uni['keywords'])})", expanded=False):
                st.write(", ".join(uni["keywords"][:24]))

    # ---- Structure + heatmap (full width, below the card) ----
    st.markdown(f"### 3D structure painted by {COND_LABEL[cond]}")

    _cpick, _ctext, _cclear = st.columns([2.2, 1.4, 0.8])
    with _cpick:
        _region_opts = ["(none)"] + [f"{s}–{e}  (aging hotspot)" for (s, e) in _aging_regions]
        _cur_idx = 0
        if focus_sel and tuple(focus_sel) in [tuple(r) for r in _aging_regions]:
            _cur_idx = [tuple(r) for r in _aging_regions].index(tuple(focus_sel)) + 1
        _pick = st.selectbox("Focus a region on the structure", _region_opts, index=_cur_idx,
                             help="Highlights the region on the 3D model and marks it on the heatmap below.")
        if _pick != "(none)":
            _pi = _region_opts.index(_pick) - 1
            _picked = tuple(_aging_regions[_pi])
            if _picked != (tuple(focus_sel) if focus_sel else None):
                st.session_state[_sel_key] = _picked
                focus_sel = _picked
    with _ctext:
        _typed = st.text_input("…or type residues", value="",
                               placeholder="e.g. 82-137 or 95",
                               help="Jump the 3D camera to any residue or range.")
        _ts = _parse_selection(_typed)
        if _ts:
            st.session_state[_sel_key] = _ts
            focus_sel = _ts
    with _cclear:
        st.write("")
        if st.button("Clear", help="Clear the 3D/heatmap selection"):
            st.session_state.pop(_sel_key, None)
            focus_sel = None

    pdb_text = get_pdb_text(uid)
    if pdb_text is None:
        st.warning(f"No AlphaFold structure available for {uid}. (Not in UP000001940 bundle or EBI.)")
    else:

        if len(sub):
            vals = sub[col_metric_col].dropna().values
            if auto_scale:
                vmax = float(np.nanpercentile(np.abs(vals), 98)) if len(vals) else 1.0
                vmax = max(vmax, 0.5)
            else:
                vmax = float(vmax_input)
        else:
            vmax = 1.0

        res2color = {int(r): diverging_color(v, vmax)
                     for r, v in zip(sub["residue"], sub[col_metric_col]) if np.isfinite(v)}

        view = py3Dmol.view(width=740, height=380)
        view.addModel(pdb_text, "pdb")
        view.setStyle({"cartoon": {"color": "#dcdcdc"}})
        for r, c in res2color.items():
            view.setStyle({"resi": str(r)}, {"cartoon": {"color": c}})

        if focus_sel:
            _fs, _fe = int(focus_sel[0]), int(focus_sel[1])
            _resi_range = f"{_fs}-{_fe}"
            view.addStyle({"resi": _resi_range},
                          {"cartoon": {"thickness": 1.0}, "stick": {"radius": 0.18, "color": "#e0a300"}})
            view.addResLabels({"resi": f"{_fs}"}, {"fontSize": 11, "backgroundColor": "#e0a300",
                                                   "backgroundOpacity": 0.7, "fontColor": "white"})
            if _fe != _fs:
                view.addResLabels({"resi": f"{_fe}"}, {"fontSize": 11, "backgroundColor": "#e0a300",
                                                      "backgroundOpacity": 0.7, "fontColor": "white"})
            view.zoomTo({"resi": _resi_range})
        else:
            view.zoomTo()
        view.setBackgroundColor("white")

        html = view._make_html()
        components.html(html, height=390, scrolling=False)
        if focus_sel:
            st.caption(f"🔎 Focused on residues **{focus_sel[0]}–{focus_sel[1]}** "
                       "(gold outline). Selection also marked on the heatmap below. Use *Clear* to reset the view.")

        fig, ax = plt.subplots(figsize=(1.9, 0.16))
        norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
        cb = mpl.colorbar.ColorbarBase(ax, cmap=LIPMS_CMAP, norm=norm, orientation="horizontal")
        cb.set_label(f"{col_metric}  (± {vmax:.2f})", fontsize=4)
        cb.ax.tick_params(labelsize=4, length=1.5, width=0.4, pad=1)
        cb.outline.set_linewidth(0.3)
        for s in ax.spines.values(): s.set_visible(False)
        st.pyplot(fig, use_container_width=False)
        plt.close(fig)
        st.caption("🔴 **Red** = positive log₂FC — region becomes more **stabilized / protected** with age · "
                   "🔵 **Blue** = negative log₂FC — region becomes more **destabilized / structurally open** · "
                   "⚪ grey = little change or not measured.")

    st.markdown("#### Residue signal — all conditions")
    if seq_len:
        cond_order = [c for c, _, _ in CONDITIONS]
        mat = np.full((len(cond_order), seq_len), np.nan, dtype=float)
        prot_all = per_res[per_res["uniprot_id"] == uid]
        for i, c in enumerate(cond_order):
            s = prot_all[prot_all["condition"] == c]
            if len(s):
                m = (s["residue"].values >= 1) & (s["residue"].values <= seq_len)
                mat[i, s["residue"].values[m] - 1] = s[col_metric_col].values[m]

        vmax_h = float(np.nanpercentile(np.abs(mat[np.isfinite(mat)]), 98)) if np.isfinite(mat).any() else 1.0
        vmax_h = max(vmax_h, 0.5)

        prot_domains = interpro[(interpro["uniprot_id"] == uid) &
                                 (interpro["domain_type"].isin(["domain", "family", "repeat"]))].copy()

        pred_df, pred_label = load_predictions()
        pred_track = None
        if len(pred_df):
            p_this = pred_df[pred_df["uniprot_id"] == uid]
            if len(p_this):
                pred_track = np.full(seq_len, np.nan, dtype=float)
                mp = (p_this["residue"].values >= 1) & (p_this["residue"].values <= seq_len)
                pred_track[p_this["residue"].values[mp] - 1] = p_this["p_destabilized"].values[mp]

        var_df = load_variants()
        var_track = None
        if len(var_df):
            v_this = var_df[var_df["worm_uniprot"] == uid]
            if len(v_this):
                var_track = np.zeros(seq_len, dtype=float)
                mv = (v_this["worm_residue"].values >= 1) & (v_this["worm_residue"].values <= seq_len)
                var_track[v_this["worm_residue"].values[mv] - 1] = v_this["n_variants_pathogenic"].values[mv]

        n_extra_rows = ((1 if pred_track is not None else 0)
                        + (1 if var_track is not None else 0)
                        + (1 if len(prot_domains) else 0))
        fig_h = 2.6 + 0.6 * n_extra_rows
        fig = plt.figure(figsize=(11, fig_h))
        h_ratios = [len(cond_order)]
        if pred_track is not None: h_ratios.append(1.2)
        if var_track is not None: h_ratios.append(1.2)
        if len(prot_domains): h_ratios.append(1.4)
        gs = fig.add_gridspec(len(h_ratios), 1, height_ratios=h_ratios, hspace=0.35)

        axh = fig.add_subplot(gs[0])
        norm = TwoSlopeNorm(vcenter=0, vmin=-vmax_h, vmax=vmax_h)
        axh.imshow(mat, aspect="auto", cmap=LIPMS_CMAP, norm=norm, interpolation="nearest",
                   extent=[0.5, seq_len+0.5, len(cond_order)-0.5, -0.5])
        axh.set_yticks(range(len(cond_order)))
        axh.set_yticklabels([COND_LABEL[c] for c in cond_order], fontsize=8)
        axh.tick_params(axis='x', labelbottom=(pred_track is None and not len(prot_domains)))
        if pred_track is None and not len(prot_domains):
            axh.set_xlabel("Residue position")
        cur = cond_order.index(cond)
        axh.axhspan(cur-0.5, cur+0.5, facecolor="none", edgecolor="#333333", linewidth=1.2)
        for s in axh.spines.values(): s.set_edgecolor("#888888")

        gs_idx = 1
        if pred_track is not None:
            axp = fig.add_subplot(gs[gs_idx], sharex=axh); gs_idx += 1
            axp.imshow(pred_track.reshape(1, -1), aspect="auto", cmap="Reds",
                       vmin=0, vmax=1, interpolation="nearest",
                       extent=[0.5, seq_len+0.5, 0.5, -0.5])
            axp.set_yticks([0]); axp.set_yticklabels([f"ML predicted\np(destabilized)\n{pred_label}"], fontsize=7)
            labels_here = (var_track is None) and (not len(prot_domains))
            axp.set_xlabel("Residue position" if labels_here else "")
            axp.tick_params(axis='x', labelbottom=labels_here)
            for s in axp.spines.values(): s.set_edgecolor("#888888")

        if var_track is not None:
            axv = fig.add_subplot(gs[gs_idx], sharex=axh); gs_idx += 1
            v_max = max(1.0, float(np.percentile(var_track[var_track > 0], 95)) if (var_track > 0).any() else 1.0)
            axv.imshow(var_track.reshape(1, -1), aspect="auto", cmap="Reds",
                        vmin=0, vmax=v_max, interpolation="nearest",
                        extent=[0.5, seq_len+0.5, 0.5, -0.5])
            axv.set_yticks([0]); axv.set_yticklabels(["Pathogenic\nvariants\n(ortholog)"], fontsize=7)
            labels_here = not len(prot_domains)
            axv.set_xlabel("Residue position" if labels_here else "")
            axv.tick_params(axis='x', labelbottom=labels_here)
            for s in axv.spines.values(): s.set_edgecolor("#888888")

        if len(prot_domains):
            axd = fig.add_subplot(gs[gs_idx], sharex=axh)
            palette = ["#8ac4d0", "#e2a25f", "#c48fbf", "#94b76c", "#d68787", "#6cb0d8", "#b0a0d0", "#c9c48c"]
            seen_col = {}
            prot_domains_sorted = prot_domains.sort_values("domain_start").reset_index(drop=True)
            LABEL_PAD_AA = max(60, int(seq_len * 0.04))
            row_ends: list[int] = []
            assignments: list[int] = []
            for _, d in prot_domains_sorted.iterrows():
                placed = False
                for i, end in enumerate(row_ends):
                    if d["domain_start"] > end + LABEL_PAD_AA:
                        row_ends[i] = int(d["domain_end"]); assignments.append(i); placed = True; break
                if not placed:
                    row_ends.append(int(d["domain_end"])); assignments.append(len(row_ends) - 1)
            n_rows_used = max(1, len(row_ends))
            for (_, d), r_idx in zip(prot_domains_sorted.iterrows(), assignments):
                key = d["interpro_id"]
                if key not in seen_col:
                    seen_col[key] = palette[len(seen_col) % len(palette)]
                axd.barh(r_idx, d["domain_end"]-d["domain_start"]+1, left=d["domain_start"]-1, height=0.55,
                         color=seen_col[key], edgecolor="#333333", linewidth=0.4)
                name = str(d["interpro_name"])[:22]
                axd.text((d["domain_start"]+d["domain_end"])/2, r_idx, name,
                         ha="center", va="center", fontsize=6.5, color="#111111")
            axd.set_ylim(n_rows_used - 0.4, -0.6); axd.set_yticks([])
            axd.set_xlim(0, seq_len)
            axd.set_xlabel("Residue")
            axd.set_title("InterPro domains", loc="left", fontsize=8)
            for s in ["top", "right", "left"]: axd.spines[s].set_visible(False)
            axd.spines["bottom"].set_edgecolor("#888888")

        if focus_sel:
            _fs, _fe = float(focus_sel[0]) - 0.5, float(focus_sel[1]) + 0.5
            for _ax in fig.axes:
                _ax.axvspan(_fs, _fe, facecolor="none", edgecolor="#e0a300", linewidth=1.6, zorder=5)

        st.pyplot(fig, use_container_width=True)
        plt.close(fig)

    # ---- Structural context sits directly below the images ----
    try:
        import structure_features as _sf
        _agg = per_res[(per_res["uniprot_id"] == uid)
                       & (per_res["condition"].isin(["day6", "day9"]))
                       & (per_res["log2fc_max_abs"].abs() > 1.0)
                       & (per_res["adj_pval_at_max_abs"] < 0.05)
                       & (per_res["adj_pval_at_max_abs"] >= 0)]
        _resis = sorted(_agg["residue"].astype(int).unique())
        _regs = []
        if _resis:
            _s = _prev = _resis[0]
            for _r in _resis[1:]:
                if _r - _prev <= 20:
                    _prev = _r
                else:
                    _regs.append((_s, _prev)); _s = _prev = _r
            _regs.append((_s, _prev))
            def _peak(rg):
                seg = _agg[(_agg["residue"] >= rg[0]) & (_agg["residue"] <= rg[1])]
                return seg["log2fc_max_abs"].abs().max()
            _regs = sorted(_regs, key=_peak, reverse=True)[:4]
        _sc = _sf.get_structural_context(uid, regions=_regs)
        if _sc.get("available"):
            st.markdown("### Structural context")
            _pp = _sc.get("per_protein", {})
            st.caption(f"AlphaFold model: {_pp.get('pct_buried','?')}% of residues buried, "
                       f"median pLDDT {_pp.get('median_plddt','?')}. "
                       f"Interprets where aging changes sit in the fold.")
            if _sc.get("regions"):
                for _r in _sc["regions"]:
                    _conf = "🟢" if _r["mean_plddt"] >= 70 else ("🟡" if _r["mean_plddt"] >= 50 else "🔴")
                    _loc = "buried core" if _r["pct_buried"] >= 50 else "surface"
                    _site = ""
                    if _r.get("near_site"):
                        _site = f" · near {_r.get('nearest_site_type','site')} ({_r.get('min_site_distance')} Å)"
                    _flag = " · ⚠️ low-confidence" if _r.get("low_confidence") else ""
                    _rc1, _rc2 = st.columns([3, 1])
                    with _rc1:
                        st.markdown(
                            f"**{_r['start']}–{_r['end']}** — {_loc}, {_r['dominant_sse']} "
                            f"· pLDDT {_r['mean_plddt']:.0f} {_conf}{_site}{_flag}"
                        )
                    with _rc2:
                        if st.button("Focus 3D", key=f"focus_btn_{uid}_{_r['start']}_{_r['end']}",
                                     help="Highlight this region on the 3D structure and heatmap"):
                            st.session_state[_sel_key] = (int(_r["start"]), int(_r["end"]))
                            st.rerun()
                    st.caption(_r["interpretation"])
            else:
                st.caption("No significant aging regions to place structurally.")
            st.caption("_AlphaFold-predicted structural properties (biotite SASA/SSE, pLDDT, "
                       "CA–CA site distance), not an experimental structure. See STRUCTURE_README.md._")
    except Exception as _e:
        pass

    # ---- Data sections above the chat ----
    _preds, _pred_label = load_predictions()
    if len(_preds):
        _p_this = _preds[_preds["uniprot_id"] == uid]
        if len(_p_this):
            st.markdown("### Machine-learning predicted aging-vulnerable residues")
            top10 = _p_this.nlargest(10, "p_destabilized")[["residue", "p_destabilized", "y_observed", "max_abs_l2fc_aging"]].copy()
            top10.columns = ["residue", "p(destabilized)", "observed?", "max |log2FC|"]
            top10["p(destabilized)"] = top10["p(destabilized)"].map(lambda v: f"{v:.3f}")
            top10["max |log2FC|"] = top10["max |log2FC|"].map(lambda v: f"{v:.2f}")
            top10["observed?"] = top10["observed?"].map({0: "no", 1: "yes"})
            st.dataframe(top10, hide_index=True, use_container_width=True, height=250)
            _perf = "OOF AUC = 0.699 · AP = 0.331" if "v2" in _pred_label else "OOF AUC = 0.679 · AP = 0.302"
            st.caption("Isotonic-calibrated LightGBM probability of |log2FC|>1 & adj-p<0.05 in aging (day6/day9). "
                       f"{_perf}. See MODEL_CARD.md for known failure modes.")

    _var_df = load_variants()
    if len(_var_df):
        _v_this = _var_df[_var_df["worm_uniprot"] == uid]
        if len(_v_this):
            st.markdown("### Disease variants (from ortholog)")
            _n_all = int(_v_this["n_variants_all"].sum())
            _n_path = int(_v_this["n_variants_pathogenic"].sum())
            _n_res_path = int((_v_this["n_variants_pathogenic"] > 0).sum())
            st.caption(f"{_n_path:,} pathogenic variants across {_n_res_path} residues "
                        f"(of {_n_all:,} total variants mapped through the top human ortholog).")
            _diseases = _v_this[_v_this["top_disease"].notna()]["top_disease"].value_counts().head(5)
            if len(_diseases):
                st.markdown("**Top diseases:**")
                for term, n in _diseases.items():
                    st.markdown(f"- {term}  <span style='color:#888;font-size:0.85em'>({n} residues)</span>", unsafe_allow_html=True)
            _top_res = _v_this.nlargest(5, "n_variants_pathogenic")
            _top_res = _top_res[_top_res["n_variants_pathogenic"] > 0]
            if len(_top_res):
                st.markdown("**Top variant-dense residues:**")
                st.dataframe(
                    _top_res[["worm_residue", "n_variants_pathogenic", "n_variants_all", "top_disease"]].rename(columns={
                        "worm_residue": "residue",
                        "n_variants_pathogenic": "path.",
                        "n_variants_all": "total",
                        "top_disease": "top disease",
                    }),
                    hide_index=True, use_container_width=True, height=180,
                )
            st.caption("_Variants mapped from human ortholog via pairwise BLOSUM62 alignment. "
                        "Homology inference, not measurements in this dataset. See VARIANTS_README.md._")

    st.markdown("---")
    st.markdown("### 🧬 Human ortholog & disease associations")
    orthologs = load_orthologs().get(uid, {})
    if not orthologs:
        st.info(f"No worm→human ortholog data available for {uid}. "
                "(Simplemine record missing; this is expected for worm-specific proteins like vitellogenins.)")
    else:
        ocol, dcol = st.columns([1, 1])
        with ocol:
            st.markdown("**Human orthologs**")
            orth = orthologs["orthologs"]
            if not orth:
                st.write("_No human orthologs identified across OMA/Compara/Panther/Inparanoid/etc._")
            else:
                orth_sorted = sorted(orth, key=lambda x: -x["n_methods"])
                for o in orth_sorted[:6]:
                    hgnc = o["hgnc_id"]
                    info = cached_hgnc(hgnc)
                    acc = info.get("accession")
                    gn  = info.get("gene_name") or "—"
                    pn  = info.get("protein_name") or ""
                    ln  = info.get("length")
                    methods_short = ", ".join(o["methods"].split(";")[:4])
                    if len(o["methods"].split(";")) > 4:
                        methods_short += ", …"
                    if acc:
                        st.markdown(
                            f"• **[{gn}](https://www.uniprot.org/uniprotkb/{acc})** — {pn} "
                            f"<span style='color:#666'>({ln} aa, [{hgnc}](https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{hgnc}))</span>",
                            unsafe_allow_html=True,
                        )
                    else:
                        st.markdown(f"• {hgnc}  — (UniProt lookup pending)")
                    st.caption(f"Supporting methods ({o['n_methods']}): {methods_short}")
                if len(orth) > 6:
                    st.caption(f"… and {len(orth)-6} more ortholog candidates.")
        with dcol:
            diseases = orthologs["diseases"]
            st.markdown(f"**Disease associations** ({len(diseases)})")
            if not diseases:
                st.write("_No disease associations in WormBase disease info._")
            else:
                uniq = {}
                for d in diseases:
                    uniq.setdefault(d["term"], d["source"])
                for term, src in list(uniq.items())[:12]:
                    st.markdown(f"• {term}  <span style='color:#888;font-size:0.85em'>({src})</span>",
                                unsafe_allow_html=True)
                if len(uniq) > 12:
                    st.caption(f"… and {len(uniq)-12} more.")
        st.caption("_Ortholog & disease associations from WormBase SimpleMine (compiled from OMA, "
                    "Ensembl Compara, Panther, Inparanoid, OrthoFinder, OrthoInspector, PhylomeDB, "
                    "Hieranoid, Roundup, TreeFam). Disease terms marked 'By Orthology' are inferred from human data._")
    st.markdown("---")
    st.markdown("## Conversation")
    try:
        import retrieval as _retr
        import agent_orchestrator as _orch
        import agents as _agsmod
        _agents_ok = True
    except Exception as _e:
        _agents_ok = False
        st.error(f"Conversation framework import failed: {_e}")

    if _agents_ok:
        st.caption("One unified chat: a coordinator routes your question to structural-biology, "
                   "aging & stress, disease-genetics and data-analyst specialists that reason over "
                   "the grounded data and return a single cited answer [S1]-[S8].")

        _roster = "".join(
            f"<span class='chip' title=\"{_agsmod.AGENTS[a]['role']}\">{_agsmod.AGENTS[a]['icon']} "
            f"{_agsmod.AGENTS[a]['name']}</span>"
            for a in _agsmod.SPECIALIST_ORDER)
        st.markdown(f"<div class='flag-row'>{_roster}</div>", unsafe_allow_html=True)

        _akey = _orch.agents_available()["has_key"]
        if _akey:
            st.markdown("<div class='status-pill live'><span class='dot'></span>"
                        "Live agent panel enabled</div>", unsafe_allow_html=True)
        else:
            st.markdown("<div class='status-pill off'><span class='dot'></span>"
                        "Offline mode - grounded data brief. Set an Anthropic API key to enable live agents."
                        "</div>", unsafe_allow_html=True)

        _AHIST = "agent_history"
        if _AHIST not in st.session_state:
            st.session_state[_AHIST] = []
        if "agent_focus_uid" not in st.session_state:
            st.session_state["agent_focus_uid"] = uid
        focus_uid = st.session_state["agent_focus_uid"]

        st.markdown("**Try a question:**")
        _suggs = [
            "What is this protein and what does it do?",
            "Which regions change with age, and where in the domain structure?",
            "Is this change aging-specific, or also seen under heat shock and polyQ?",
            "Is the aging region buried in the fold, and could it affect stability?",
            "Are there human disease variants near the aging region?",
            "Which residues does the ML model prioritize?",
        ]
        _prefill = None
        _sc = st.columns(2)
        for _i, _s in enumerate(_suggs):
            if _sc[_i % 2].button(_s, key=f"convsug_{_i}", use_container_width=True):
                _prefill = _s

        for _m in st.session_state[_AHIST]:
            with st.chat_message(_m["role"]):
                if _m.get("trace"):
                    _chips = " ".join(f"{t['icon']} {t['name']}" for t in _m["trace"])
                    st.caption(f"Agents consulted: {_chips}")
                st.markdown(_m["content"])

        _ain = st.chat_input("Ask about this protein or the aging dataset...") or _prefill

        if _ain:
            _new_uid = _resolve_uid_from_text(_ain, fallback=focus_uid)
            if _new_uid and _new_uid != focus_uid:
                focus_uid = _new_uid
                st.session_state["agent_focus_uid"] = _new_uid

            st.session_state[_AHIST].append({"role": "user", "content": _ain})
            with st.chat_message("user"):
                st.markdown(_ain)

            _actx = _retr.get_protein_context(focus_uid)
            with st.chat_message("assistant"):
                _route_ph = st.empty()
                _status_ph = st.empty()
                _ans_ph = st.empty()
                _acc = ""
                _trace = []
                _final = None
                _done = []
                try:
                    for _ev in _orch.ask_agents(_ain, _actx, uid=focus_uid):
                        if isinstance(_ev, str):
                            _acc += _ev
                            _ans_ph.markdown(_acc + "\u258c")
                        elif _ev.get("event") == "route":
                            _names = " . ".join(
                                f"{_agsmod.AGENTS[a]['icon']} {_agsmod.AGENTS[a]['name']}"
                                for a in _ev["specialists"])
                            _route_ph.caption(f"Routing to: {_names}")
                        elif _ev.get("event") == "specialist_done":
                            _done.append(f"{_ev['icon']} {_ev['name']} {'ok' if _ev['ok'] else 'warn'}")
                            _status_ph.caption("  ".join(_done))
                        elif _ev.get("event") == "final":
                            _final = _ev
                    _ans_ph.markdown(_acc)
                    _status_ph.empty()
                except Exception as _e:
                    _acc = f"Error: {_e}"
                    _ans_ph.markdown(_acc)

                if _final:
                    _trace = _final.get("trace", [])
                    if _final.get("usage"):
                        _u = _final["usage"]
                        st.caption(f"Model: {_final.get('model','?')} - synthesis "
                                   f"in={_u.get('input_tokens','?')} out={_u.get('output_tokens','?')} tokens - "
                                   f"{len(_trace)} specialists")
                    if _trace:
                        with st.expander("Agent contributions (what each specialist said)"):
                            for _t in _trace:
                                _src = ", ".join(_t.get("sources", []))
                                st.markdown(f"**{_t['icon']} {_t['name']}**  -  sources {_src}")
                                if _t.get("error"):
                                    st.caption(f"warn {_t['error']}")
                                elif _t.get("excerpt"):
                                    st.caption(_t["excerpt"] + ("..." if _t.get("chars", 0) > 600 else ""))
                                st.markdown("---")

            st.session_state[_AHIST].append(
                {"role": "assistant", "content": _acc,
                 "trace": [{"icon": t["icon"], "name": t["name"]} for t in _trace]})

        if focus_uid != uid:
            _fp = prot_feat[prot_feat["uniprot_id"] == focus_uid]
            _fgene = _fp.iloc[0]["gene_symbol"] if len(_fp) and pd.notna(_fp.iloc[0]["gene_symbol"]) else focus_uid
            if st.button(f"Load {_fgene} ({focus_uid}) in the viewer", use_container_width=True):
                st.session_state["protein_query"] = focus_uid
                st.rerun()
        if st.session_state.get(_AHIST):
            if st.button("Clear conversation", use_container_width=True):
                st.session_state[_AHIST] = []
                st.rerun()


    # ---------------- Footer ----------------
    st.markdown("---")
    st.caption("Dataset: 6,823 proteins × 9 conditions from TMT-LiP-MS in C. elegans. "
               "Structures from AlphaFold UP000001940 (94.7% coverage). "
               "Every claim in chat is cited to structured sources [S1]–[S8].")
