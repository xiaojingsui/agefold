"""
Discovery — per-protein composite score for ranking worm proteins by combined
aging + disease-variant signal.

The score is deliberately a HEURISTIC ranker, not a p-value. Weights are
tunable so the Streamlit tab can let the user explore what "interesting"
means to them.

Components (all per-protein, normalized per residue where noted):
  aging_hotspot_max      max |log2FC| in aging (day6/day9) for significant peptides
  aging_load_per_aa      sig aging peptides × mean |log2FC|, divided by protein length
  ml_predicted_high_frac fraction of residues with predicted p(destab) > 0.5
  variant_load_per_aa    # pathogenic variants mapped from human ortholog / protein length
  aging_variant_overlap  # pathogenic variants that fall inside any aging top-region
  alignment_confidence   0..1 quality of worm↔human alignment (identity clipped at 0.5, then rescaled)

composite = alignment_confidence × (
    w_hot * aging_hotspot_max
  + w_load * aging_load_per_aa
  + w_ml * ml_predicted_high_frac
  + w_var * variant_load_per_aa
  + w_ovl * aging_variant_overlap
)

Defaults favor overlap (the actual novel signal) over any single component alone.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd

APP_DIR = Path(__file__).parent
DATA_DIR = APP_DIR / "data"

# Default weights — tuned so top hits look reasonable (actin/myosin/HSPs high,
# constitutively-expressed housekeeping middle, no-signal proteins low)
DEFAULT_WEIGHTS = {
    "w_hot": 1.0,     # peak intensity — captures "at least one strong hotspot"
    "w_load": 5.0,    # aging load per aa — a whole-protein hit not just one spike
    "w_ml": 2.0,      # ML confidence — corroborates model believes the hits
    "w_var": 0.5,     # variant density per aa — human disease relevance
    "w_ovl": 3.0,     # overlap bonus — the actually novel signal
}

AGING_CONDITIONS = ("day6", "day9")
SIG_LOG2FC = 1.0
SIG_ADJPVAL = 0.05


def _clip_alignment_confidence(identity: float | None) -> float:
    """Alignment identity < 0.3 → 0 (unusable); 0.3-0.6 → linearly rescaled to 0-1; >=0.6 → 1."""
    if identity is None or pd.isna(identity):
        return 0.0
    if identity < 0.3:
        return 0.0
    if identity >= 0.6:
        return 1.0
    return (identity - 0.3) / 0.3


def compute_scores(weights: dict[str, float] | None = None) -> pd.DataFrame:
    """Compute per-protein composite scores and all component metrics.
    Returns one row per measured protein, sorted by composite descending."""
    W = {**DEFAULT_WEIGHTS, **(weights or {})}

    prot_feat = pd.read_parquet(DATA_DIR / "protein_features.parquet").drop_duplicates("uniprot_id")
    peptides = pd.read_parquet(DATA_DIR / "lipms_peptide_unified.parquet")
    per_res = pd.read_parquet(DATA_DIR / "lipms_per_residue_agg.parquet")

    # Restrict scoring to proteins that actually have LiP-MS measurements.
    # Unmeasured proteins in protein_features cannot have aging signal, so they
    # would sit at composite=0 and only inflate the "of N total" denominator.
    measured_uids = set(per_res["uniprot_id"].unique())
    prot_feat = prot_feat[prot_feat["uniprot_id"].isin(measured_uids)].copy()

    # ML predictions (may be v1 or v2)
    preds = None
    for name in ("residue_predictions_v2.parquet", "residue_predictions.parquet"):
        p = DATA_DIR / name
        if p.exists():
            preds = pd.read_parquet(p); break

    # Variants
    var_p = DATA_DIR / "variant_density.parquet"
    variants = pd.read_parquet(var_p) if var_p.exists() else pd.DataFrame(
        columns=["worm_uniprot", "worm_residue", "n_variants_pathogenic", "n_variants_all", "top_disease", "alignment_identity"])

    # Filter peptides to significant aging entries
    aging_pep = peptides[
        peptides["condition"].isin(AGING_CONDITIONS)
        & (peptides["adj_pval"] < SIG_ADJPVAL)
        & (peptides["adj_pval"] >= 0)
        & (peptides["log2fc"].abs() > SIG_LOG2FC)
    ].copy()

    # aging_hotspot_max: max |log2FC| per protein in significant aging peptides
    if len(aging_pep):
        hotspot = aging_pep.assign(abs_l2=aging_pep["log2fc"].abs()).groupby("uniprot_id")["abs_l2"].max()
    else:
        hotspot = pd.Series(dtype=float)

    # aging_load: n_sig × mean(|log2FC|); we'll normalize by seq_length below
    if len(aging_pep):
        agg = aging_pep.assign(abs_l2=aging_pep["log2fc"].abs()).groupby("uniprot_id").agg(
            n_sig=("peptide_sequence", "nunique"),
            mean_abs_l2=("abs_l2", "mean"),
        )
        agg["aging_load"] = agg["n_sig"] * agg["mean_abs_l2"]
    else:
        agg = pd.DataFrame(columns=["n_sig", "mean_abs_l2", "aging_load"])

    # ml_predicted_high_frac: fraction of residues with p>0.5
    if preds is not None and len(preds):
        ml_frac = preds.groupby("uniprot_id").apply(lambda g: float((g["p_destabilized"] > 0.5).mean()))
        ml_frac.name = "ml_predicted_high_frac"
    else:
        ml_frac = pd.Series(dtype=float, name="ml_predicted_high_frac")

    # variant_load: total pathogenic per protein, plus alignment identity
    if len(variants):
        vagg = variants.groupby("worm_uniprot").agg(
            variant_load=("n_variants_pathogenic", "sum"),
            n_res_pathogenic=("n_variants_pathogenic", lambda s: int((s > 0).sum())),
            alignment_identity=("alignment_identity", "max"),
        ).rename_axis("uniprot_id")
    else:
        vagg = pd.DataFrame(columns=["variant_load", "n_res_pathogenic", "alignment_identity"])

    # Aging regions (residue-level) — significant aging peptides mapped to residue positions
    # Then compute aging_variant_overlap = # pathogenic variants within aging residues
    aging_residues = per_res[
        (per_res["condition"].isin(AGING_CONDITIONS))
        & (per_res["log2fc_max_abs"].abs() > SIG_LOG2FC)
        & (per_res["adj_pval_at_max_abs"] < SIG_ADJPVAL)
        & (per_res["adj_pval_at_max_abs"] >= 0)
    ][["uniprot_id", "residue"]].drop_duplicates()

    if len(variants) and len(aging_residues):
        merged = variants.merge(aging_residues,
                                 left_on=["worm_uniprot", "worm_residue"],
                                 right_on=["uniprot_id", "residue"], how="inner")
        overlap = merged.groupby("worm_uniprot")["n_variants_pathogenic"].sum().rename_axis("uniprot_id")
        overlap.name = "aging_variant_overlap"
    else:
        overlap = pd.Series(dtype=float, name="aging_variant_overlap")

    # Assemble
    df = prot_feat[["uniprot_id", "gene_symbol", "seq_length"]].copy()
    df["aging_hotspot_max"] = df["uniprot_id"].map(hotspot).fillna(0.0)
    df["n_sig_aging_peptides"] = df["uniprot_id"].map(agg["n_sig"]).fillna(0).astype(int) if len(agg) else 0
    df["aging_load"] = df["uniprot_id"].map(agg["aging_load"]).fillna(0.0) if len(agg) else 0.0
    df["aging_load_per_aa"] = df["aging_load"] / df["seq_length"].replace(0, np.nan)
    df["aging_load_per_aa"] = df["aging_load_per_aa"].fillna(0.0)
    df["ml_predicted_high_frac"] = df["uniprot_id"].map(ml_frac).fillna(0.0)
    df["variant_load"] = df["uniprot_id"].map(vagg["variant_load"] if len(vagg) else pd.Series(dtype=float)).fillna(0).astype(int)
    df["variant_load_per_aa"] = df["variant_load"] / df["seq_length"].replace(0, np.nan)
    df["variant_load_per_aa"] = df["variant_load_per_aa"].fillna(0.0)
    df["n_res_pathogenic"] = df["uniprot_id"].map(vagg["n_res_pathogenic"] if len(vagg) else pd.Series(dtype=int)).fillna(0).astype(int)
    df["alignment_identity"] = df["uniprot_id"].map(vagg["alignment_identity"] if len(vagg) else pd.Series(dtype=float))
    df["aging_variant_overlap"] = df["uniprot_id"].map(overlap).fillna(0).astype(int)
    df["alignment_confidence"] = df["alignment_identity"].apply(_clip_alignment_confidence)

    # Top disease per protein (most common top_disease across variant-carrying residues)
    if len(variants):
        td = variants[variants["top_disease"].notna() & (variants["n_variants_pathogenic"] > 0)]
        if len(td):
            top_dis = td.groupby("worm_uniprot")["top_disease"].agg(
                lambda s: s.value_counts().index[0] if len(s) else None).rename_axis("uniprot_id")
            df["top_disease"] = df["uniprot_id"].map(top_dis)
        else:
            df["top_disease"] = None
    else:
        df["top_disease"] = None

    # Composite (alignment confidence multiplies the whole variant-relevant portion)
    variant_component = (
        W["w_var"] * df["variant_load_per_aa"]
        + W["w_ovl"] * df["aging_variant_overlap"]
    )
    aging_component = (
        W["w_hot"] * df["aging_hotspot_max"]
        + W["w_load"] * df["aging_load_per_aa"]
        + W["w_ml"] * df["ml_predicted_high_frac"]
    )
    df["composite"] = aging_component + df["alignment_confidence"] * variant_component

    df = df.sort_values("composite", ascending=False).reset_index(drop=True)
    return df


if __name__ == "__main__":
    df = compute_scores()
    print(f"Scored {len(df):,} proteins")
    print(f"Top 20:")
    cols = ["uniprot_id", "gene_symbol", "seq_length", "composite",
            "aging_hotspot_max", "aging_load_per_aa", "ml_predicted_high_frac",
            "variant_load", "aging_variant_overlap", "alignment_identity", "top_disease"]
    print(df[cols].head(20).to_string(index=False))
