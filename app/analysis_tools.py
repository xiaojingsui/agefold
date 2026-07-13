"""
analysis_tools.py — dataset-level analysis functions for the Data Analyst agent.

Bounded, precomputed-table-backed queries that answer cross-protein and
cross-condition questions (not just single-protein). Each returns a compact dict
with a `summary` string the agent can quote, plus structured `data` for the UI.

All numbers come from the same parquets the viewer uses:
  lipms_per_residue_agg.parquet, discovery_scores.parquet, variant_density.parquet,
  structural_features.parquet, site_proximity.parquet.
"""
from __future__ import annotations
import os
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd

_DATA = Path(__file__).parent / "data"
_cache: dict[str, Any] = {}

# Condition → group, mirrors streamlit_app.CONDITIONS
COND_GROUP = {
    "day6": "aging", "day9": "aging",
    "hs": "stress",
    "q35": "polyQ", "q40": "polyQ",
    "myosin_ts_15": "ts-mutant", "myosin_ts_25": "ts-mutant",
    "paramyosin_ts_15": "ts-mutant", "paramyosin_ts_25": "ts-mutant",
}
COND_LABEL = {
    "day6": "WT day 6", "day9": "WT day 9", "hs": "Heat shock",
    "q35": "Q35 (polyQ)", "q40": "Q40 (polyQ)",
    "myosin_ts_15": "myosin-ts 15°C", "myosin_ts_25": "myosin-ts 25°C",
    "paramyosin_ts_15": "paramyosin-ts 15°C", "paramyosin_ts_25": "paramyosin-ts 25°C",
}
SIG_LOG2FC = 1.0
SIG_PVAL = 0.05


def _load(name: str) -> pd.DataFrame:
    if name not in _cache:
        p = _DATA / name
        _cache[name] = pd.read_parquet(p) if p.exists() else pd.DataFrame()
    return _cache[name]


def _sig(df: pd.DataFrame) -> pd.DataFrame:
    """Significant rows: |log2FC|>1 AND 0<=adj_p<0.05."""
    return df[(df["log2fc_max_abs"].abs() > SIG_LOG2FC)
              & (df["adj_pval_at_max_abs"] < SIG_PVAL)
              & (df["adj_pval_at_max_abs"] >= 0)]


# ---------------------------------------------------------------------------
# 1. Per-protein cross-condition comparison
# ---------------------------------------------------------------------------
def compare_conditions(uid: str) -> dict[str, Any]:
    """Peak |log2FC| and significant-peptide count per condition for one protein,
    grouped into aging / stress / polyQ / ts-mutant. Answers 'is this
    aging-specific or a generic stress response?'"""
    per = _load("lipms_per_residue_agg.parquet")
    d = per[per["uniprot_id"] == uid]
    if not len(d):
        return {"available": False, "uniprot_id": uid, "summary": f"No LiP-MS data for {uid}."}
    rows = []
    for cond, sub in d.groupby("condition", observed=True):
        sig = _sig(sub)
        rows.append({
            "condition": cond, "label": COND_LABEL.get(cond, cond),
            "group": COND_GROUP.get(cond, "?"),
            "peak_abs_log2fc": round(float(sub["log2fc_max_abs"].abs().max()), 3),
            "n_sig_peptides": int(sig["n_peptides"].sum()) if len(sig) else 0,
            "n_sig_residues": int(sig["residue"].nunique()) if len(sig) else 0,
        })
    tab = pd.DataFrame(rows)
    # group-level peaks
    grp = tab.groupby("group")["peak_abs_log2fc"].max().to_dict()
    aging_peak = grp.get("aging", 0.0)
    other_peak = max([v for g, v in grp.items() if g != "aging"], default=0.0)
    if aging_peak >= 1.0 and aging_peak >= 1.5 * max(other_peak, 1e-6):
        verdict = "aging-specific (aging peak clearly exceeds stress/polyQ/ts)"
    elif aging_peak >= 1.0 and other_peak >= 1.0:
        verdict = "shared — the region also responds to non-aging stress (generic destabilization)"
    elif aging_peak < 1.0:
        verdict = "little aging signal in this protein"
    else:
        verdict = "mixed"
    summary = (f"{uid}: aging peak |log2FC|={aging_peak:.2f}, "
               f"strongest non-aging peak={other_peak:.2f} → {verdict}. "
               f"Group peaks: " + ", ".join(f"{g} {v:.2f}" for g, v in sorted(grp.items())))
    return {"available": True, "uniprot_id": uid, "group_peaks": grp,
            "verdict": verdict, "table": tab.to_dict("records"), "summary": summary}


# ---------------------------------------------------------------------------
# 2. Top movers in a condition
# ---------------------------------------------------------------------------
def top_movers(condition: str = "day6", n: int = 15) -> dict[str, Any]:
    """Proteins with the largest peak |log2FC| in a given condition."""
    per = _load("lipms_per_residue_agg.parquet")
    d = per[per["condition"] == condition]
    if not len(d):
        return {"available": False, "condition": condition,
                "summary": f"No data for condition '{condition}'. Valid: {sorted(COND_LABEL)}."}
    g = (d.groupby(["uniprot_id", "gene_symbol"], observed=True)["log2fc_max_abs"]
         .apply(lambda x: x.abs().max()).reset_index(name="peak_abs_log2fc")
         .sort_values("peak_abs_log2fc", ascending=False).head(n))
    recs = g.to_dict("records")
    top = ", ".join(f"{r['gene_symbol']} ({r['peak_abs_log2fc']:.2f})" for r in recs[:5])
    summary = (f"Top movers in {COND_LABEL.get(condition, condition)}: {top}"
               + (f" … (top {n} returned)" if n > 5 else ""))
    return {"available": True, "condition": condition, "label": COND_LABEL.get(condition, condition),
            "table": recs, "summary": summary}


# ---------------------------------------------------------------------------
# 3. Aging-specific hits across the proteome
# ---------------------------------------------------------------------------
def aging_specific_hits(n: int = 20, min_ratio: float = 1.5) -> dict[str, Any]:
    """Proteins whose aging peak |log2FC| is high AND clearly exceeds their peak
    under stress/polyQ/ts (aging-specific, not generic destabilization)."""
    per = _load("lipms_per_residue_agg.parquet")
    if not len(per):
        return {"available": False, "summary": "per-residue table unavailable."}
    per = per.copy()
    per["group"] = per["condition"].map(COND_GROUP)
    peak = (per.groupby(["uniprot_id", "gene_symbol", "group"], observed=True)["log2fc_max_abs"]
            .apply(lambda x: x.abs().max()).reset_index(name="peak"))
    wide = peak.pivot_table(index=["uniprot_id", "gene_symbol"], columns="group",
                            values="peak", aggfunc="max").fillna(0.0)
    for g in ["aging", "stress", "polyQ", "ts-mutant"]:
        if g not in wide.columns:
            wide[g] = 0.0
    wide["other_max"] = wide[["stress", "polyQ", "ts-mutant"]].max(axis=1)
    wide["aging_specificity"] = wide["aging"] / wide["other_max"].clip(lower=0.3)
    hits = wide[(wide["aging"] >= 1.0) & (wide["aging"] >= min_ratio * wide["other_max"].clip(lower=0.3))]
    hits = hits.sort_values("aging", ascending=False).head(n).reset_index()
    recs = [{"uniprot_id": r["uniprot_id"], "gene_symbol": r["gene_symbol"],
             "aging_peak": round(r["aging"], 2), "other_peak": round(r["other_max"], 2),
             "specificity": round(r["aging_specificity"], 2)} for _, r in hits.iterrows()]
    top = ", ".join(f"{r['gene_symbol']} (aging {r['aging_peak']} vs other {r['other_peak']})" for r in recs[:6])
    summary = (f"{len(recs)} strongest aging-specific proteins (aging peak ≥1.0 and ≥{min_ratio}× "
               f"the stress/polyQ/ts peak): {top}")
    return {"available": True, "n": len(recs), "min_ratio": min_ratio,
            "table": recs, "summary": summary}


# ---------------------------------------------------------------------------
# 4. Structural enrichment summary (whole dataset)
# ---------------------------------------------------------------------------
def structural_enrichment_summary() -> dict[str, Any]:
    """Are aging hotspots enriched in buried/ordered/near-site regions vs the
    measured background? (the whole-dataset structural finding)."""
    sf = _load("structural_features.parquet")
    per = _load("lipms_per_residue_agg.parquet")
    prox = _load("site_proximity.parquet")
    if not len(sf) or not len(per):
        return {"available": False, "summary": "structural/per-residue tables unavailable."}
    aging = _sig(per[per["condition"].isin(["day6", "day9"])])[["uniprot_id", "residue"]].drop_duplicates()
    h = aging.merge(sf, on=["uniprot_id", "residue"], how="inner")
    bg = sf.merge(per[["uniprot_id", "residue"]].drop_duplicates(), on=["uniprot_id", "residue"], how="inner")
    out = {
        "available": True,
        "n_hotspot_residues": int(len(h)),
        "hotspot_pct_buried": round(100 * h["buried"].mean(), 1),
        "background_pct_buried": round(100 * bg["buried"].mean(), 1),
        "hotspot_median_plddt": round(float(h["plddt"].median()), 1),
        "background_median_plddt": round(float(bg["plddt"].median()), 1),
    }
    if len(prox):
        hs = aging.merge(prox, on=["uniprot_id", "residue"], how="inner")
        bgs = per[["uniprot_id", "residue"]].drop_duplicates().merge(prox, on=["uniprot_id", "residue"], how="inner")
        out["hotspot_near_site_pct"] = round(100 * hs["near_site"].mean(), 1)
        out["background_near_site_pct"] = round(100 * bgs["near_site"].mean(), 1)
        out["near_site_enrichment"] = round(out["hotspot_near_site_pct"] / max(out["background_near_site_pct"], 1e-6), 2)
    out["summary"] = (
        f"Across {out['n_hotspot_residues']:,} aging-hotspot residues: {out['hotspot_pct_buried']}% buried "
        f"vs {out['background_pct_buried']}% background (similar); "
        + (f"but {out.get('near_site_enrichment','?')}× enriched near functional sites "
           f"({out.get('hotspot_near_site_pct')}% vs {out.get('background_near_site_pct')}%)."
           if "near_site_enrichment" in out else "")
    )
    return out


# ---------------------------------------------------------------------------
# 5. Discovery top hits (disease-linked candidates)
# ---------------------------------------------------------------------------
def discovery_top(n: int = 10, require_disease: bool = True) -> dict[str, Any]:
    """Top composite-scored candidates from the discovery dashboard."""
    dsc = _load("discovery_scores.parquet")
    if not len(dsc):
        return {"available": False, "summary": "discovery_scores table unavailable."}
    d = dsc.copy()
    if require_disease:
        d = d[d["top_disease"].notna() & (d["top_disease"] != "")]
    d = d.sort_values("composite", ascending=False).head(n)
    recs = [{"uniprot_id": r["uniprot_id"], "gene_symbol": r["gene_symbol"],
             "composite": round(float(r["composite"]), 1),
             "top_disease": r["top_disease"]} for _, r in d.iterrows()]
    top = "; ".join(f"{r['gene_symbol']} → {r['top_disease']} ({r['composite']})" for r in recs[:5])
    summary = f"Top discovery candidates by composite score: {top}"
    return {"available": True, "table": recs, "summary": summary}


# Registry the orchestrator can expose to the Data Analyst.
DATASET_TOOLS = {
    "compare_conditions": compare_conditions,
    "top_movers": top_movers,
    "aging_specific_hits": aging_specific_hits,
    "structural_enrichment_summary": structural_enrichment_summary,
    "discovery_top": discovery_top,
}
