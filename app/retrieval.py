"""
Retrieval layer for the protein-page chat.

`get_protein_context(uid)` assembles everything the app knows about a protein
into a structured dict, ready to be templated into an LLM prompt or rendered
as an offline report card. Every retrieved fact is tagged with a source id
(S1, S2, ...) so downstream prompts can enforce citation.
"""
from __future__ import annotations
from pathlib import Path
from typing import Any
import json
import pandas as pd
import numpy as np

APP_DIR = Path(__file__).parent
DATA_DIR = APP_DIR / "data"

# Lazy per-process cache of the parquet layer
_cache: dict[str, Any] = {}


def _p(name: str) -> pd.DataFrame:
    """Cached parquet reader."""
    if name not in _cache:
        _cache[name] = pd.read_parquet(DATA_DIR / name)
    return _cache[name]


def _safe(val, default="—"):
    if val is None: return default
    try:
        if isinstance(val, float) and np.isnan(val): return default
    except Exception: pass
    if isinstance(val, (list, dict)) and len(val) == 0: return default
    return val


def _peptide_summary(uid: str, peptides: pd.DataFrame) -> dict[str, Any]:
    """Per-condition peptide counts and top-|log2FC| residues."""
    slc = peptides[peptides["uniprot_id"] == uid]
    if len(slc) == 0:
        return {"per_condition": {}, "n_peptides_total": 0}
    per_condition = {}
    for cond, g in slc.groupby("condition"):
        top = g.reindex(g["log2fc"].abs().sort_values(ascending=False).index).head(3)
        per_condition[cond] = {
            "n_peptides": int(len(g)),
            "n_significant": int(((g["log2fc"].abs() > 1) & (g["adj_pval"] < 0.05)).sum()),
            "top_peptides": [
                {
                    "sequence": r["peptide_sequence"],
                    "start": int(r["start_pos"]) if pd.notna(r["start_pos"]) else None,
                    "end": int(r["end_pos"]) if pd.notna(r["end_pos"]) else None,
                    "log2fc": round(float(r["log2fc"]), 3),
                    "adj_pval": round(float(r["adj_pval"]), 4) if pd.notna(r["adj_pval"]) else None,
                    "type": r["peptide_type"],
                }
                for _, r in top.iterrows()
            ],
        }
    return {"per_condition": per_condition, "n_peptides_total": int(slc["peptide_sequence"].nunique())}


def _lipms_aging_signal(uid: str, peptides: pd.DataFrame) -> dict[str, Any]:
    """Focused view: which regions of the protein change with age."""
    aging = peptides[(peptides["uniprot_id"] == uid) &
                     (peptides["condition"].isin(["day6", "day9"])) &
                     (peptides["adj_pval"] < 0.05) &
                     (peptides["adj_pval"] >= 0) &
                     (peptides["log2fc"].abs() > 1)].copy()
    if len(aging) == 0:
        return {"n_significant_aging_peptides": 0, "top_regions": [], "conditions": {}}

    # Cluster into contiguous regions (merge peptides within 20 aa)
    aging["mid"] = ((aging["start_pos"].fillna(0) + aging["end_pos"].fillna(0)) / 2).astype(int)
    aging = aging.sort_values("mid").reset_index(drop=True)
    regions: list[dict] = []
    for _, r in aging.iterrows():
        st, en = int(r["start_pos"]), int(r["end_pos"])
        if regions and st <= regions[-1]["end"] + 20:
            regions[-1]["start"] = min(regions[-1]["start"], st)
            regions[-1]["end"] = max(regions[-1]["end"], en)
            regions[-1]["max_abs_log2fc"] = max(regions[-1]["max_abs_log2fc"], abs(float(r["log2fc"])))
            regions[-1]["n_peptides"] += 1
        else:
            regions.append({
                "start": st, "end": en, "n_peptides": 1,
                "max_abs_log2fc": abs(float(r["log2fc"])),
            })
    # Top 5 regions by max |log2FC|
    regions = sorted(regions, key=lambda r: -r["max_abs_log2fc"])[:5]
    for r in regions:
        r["max_abs_log2fc"] = round(r["max_abs_log2fc"], 3)

    cond_counts = aging.groupby("condition")["peptide_sequence"].nunique().to_dict()
    return {
        "n_significant_aging_peptides": int(len(aging)),
        "top_regions": regions,
        "conditions": {c: int(v) for c, v in cond_counts.items()},
    }


def _annotate_regions_with_domains(regions: list[dict], domains: list[dict]) -> list[dict]:
    """For each region, list the domains it overlaps (by ≥50% or ≥30 aa overlap).
    Attaches an 'in_domains' key so the LLM doesn't have to re-derive interval math."""
    for r in regions:
        r_start, r_end = int(r["start"]), int(r["end"])
        r_len = r_end - r_start + 1
        overlaps = []
        for d in domains:
            d_start, d_end = int(d["start"]), int(d["end"])
            ov = max(0, min(r_end, d_end) - max(r_start, d_start) + 1)
            if ov <= 0: continue
            frac = ov / r_len
            if frac >= 0.5 or ov >= 30:
                overlaps.append({
                    "interpro_id": d["interpro_id"],
                    "name": d["name"],
                    "type": d["type"],
                    "span": f"{d_start}-{d_end}",
                    "overlap_aa": int(ov),
                    "fraction_of_region": round(frac, 2),
                })
        # Sort by overlap size descending, keep top 2
        r["in_domains"] = sorted(overlaps, key=lambda x: -x["overlap_aa"])[:2]
    return regions


def _ml_predictions(uid: str) -> dict[str, Any]:
    """Top predicted aging-vulnerable residues from the ML model."""
    for name, ver in [("residue_predictions_v2.parquet", "v2"),
                       ("residue_predictions.parquet", "v1")]:
        path = DATA_DIR / name
        if not path.exists(): continue
        df = pd.read_parquet(path, filters=[("uniprot_id", "==", uid)])
        if len(df) == 0: continue
        top = df.nlargest(10, "p_destabilized")
        return {
            "model_version": ver,
            "n_residues_scored": int(len(df)),
            "n_high_confidence": int((df["p_destabilized"] > 0.5).sum()),
            "top_residues": [
                {
                    "residue": int(r["residue"]),
                    "p_destabilized": round(float(r["p_destabilized"]), 3),
                    "observed": bool(int(r["y_observed"])),
                    "max_abs_log2fc_observed": round(float(r["max_abs_l2fc_aging"]), 3),
                }
                for _, r in top.iterrows()
            ],
        }
    return {"model_version": None, "n_residues_scored": 0, "n_high_confidence": 0, "top_residues": []}


def _uniprot_enrichment(uid: str) -> dict[str, Any]:
    """Read the enrichment SQLite cache (populated by the app during normal browsing)."""
    import sqlite3
    path = DATA_DIR / "enrichment_cache.sqlite"
    if not path.exists():
        return {}
    conn = sqlite3.connect(str(path))
    try:
        row = conn.execute("SELECT payload FROM enrichment WHERE uniprot_id=? AND source='uniprot'",
                            (uid,)).fetchone()
        if not row: return {}
        payload = json.loads(row[0])
        # Trim to essentials
        return {
            "recommended_name": payload.get("recommended_name"),
            "function_text": payload.get("function") or payload.get("function_text"),
            "keywords": (payload.get("keywords") or [])[:8],
            "subcellular_location": (payload.get("subcellular_location") or [])[:4],
            "go_bp": [f"{g['id']} {g['name']}" for g in (payload.get("go_bp") or [])[:5]],
            "go_mf": [f"{g['id']} {g['name']}" for g in (payload.get("go_mf") or [])[:5]],
            "go_cc": [f"{g['id']} {g['name']}" for g in (payload.get("go_cc") or [])[:3]],
        }
    finally:
        conn.close()


def _domains(uid: str) -> list[dict[str, Any]]:
    ipr = _p("interpro_domains.parquet")
    d = ipr[ipr["uniprot_id"] == uid].copy()
    if len(d) == 0: return []
    d["size"] = d["domain_end"].astype(int) - d["domain_start"].astype(int) + 1
    d = d.sort_values("size", ascending=False).head(10)
    return [
        {
            "interpro_id": r["interpro_id"],
            "name": r["interpro_name"],
            "type": r["domain_type"],
            "start": int(r["domain_start"]),
            "end": int(r["domain_end"]),
        }
        for _, r in d.iterrows()
    ]


def _orthologs(uid: str) -> dict[str, Any]:
    df = _p("orthologs_and_disease.parquet")
    row = df[df["worm_uniprot"] == uid]
    if len(row) == 0:
        return {"orthologs": [], "diseases": [], "concise_description": None}
    r = row.iloc[0]
    try:
        orth = json.loads(r["orthologs"]) if r["orthologs"] else []
    except Exception:
        orth = []
    try:
        dis = json.loads(r["diseases"]) if r["diseases"] else []
    except Exception:
        dis = []
    # Top 3 orthologs by method count (n_methods is a pre-stored int; methods is a semicolon-joined string)
    orth = sorted(orth, key=lambda o: -int(o.get("n_methods", 0)))[:3]
    # Deduplicate disease terms
    seen = set()
    unique_dis = []
    for d in dis:
        term = d.get("term")
        if term and term not in seen:
            seen.add(term); unique_dis.append(d)
    unique_dis = unique_dis[:10]
    return {
        "orthologs": [{"hgnc_id": o["hgnc_id"], "n_methods": int(o.get("n_methods", 0))} for o in orth],
        "diseases": [d["term"] for d in unique_dis],
        "concise_description": r.get("concise_description"),
    }


def _variant_overlap(uid: str, aging_regions: list[dict]) -> dict[str, Any]:
    """For each aging region, count pathogenic variants that fall within it
    in the mapped ortholog. Returns per-region overlap plus a top-diseases summary."""
    path = DATA_DIR / "variant_density.parquet"
    if not path.exists():
        return {"available": False, "regions": []}
    vdf = pd.read_parquet(path, filters=[("worm_uniprot", "==", uid)])
    if len(vdf) == 0:
        return {"available": True, "n_pathogenic_total": 0, "n_residues_pathogenic": 0,
                "top_diseases": [], "regions": []}
    per_region = []
    for r in aging_regions:
        st, en = int(r["start"]), int(r["end"])
        rslice = vdf[(vdf["worm_residue"] >= st) & (vdf["worm_residue"] <= en)]
        n_path = int(rslice["n_variants_pathogenic"].sum())
        n_all = int(rslice["n_variants_all"].sum())
        top_dis = [d for d in rslice["top_disease"].dropna().unique() if d][:3]
        per_region.append({"start": st, "end": en, "n_pathogenic_in_region": n_path,
                            "n_variants_all_in_region": n_all, "top_diseases_in_region": top_dis})
    # Overall top diseases across the whole protein
    all_dis = [d for d in vdf["top_disease"].dropna() if d]
    from collections import Counter
    top_disease_overall = [t for t, _ in Counter(all_dis).most_common(5)]
    return {
        "available": True,
        "n_pathogenic_total": int(vdf["n_variants_pathogenic"].sum()),
        "n_residues_pathogenic": int((vdf["n_variants_pathogenic"] > 0).sum()),
        "n_variants_all_total": int(vdf["n_variants_all"].sum()),
        "top_diseases": top_disease_overall,
        "regions": per_region,
    }


def _structural_context(uid: str, aging_regions: list) -> dict[str, Any]:
    """AlphaFold-derived structural context for the aging regions (buried/exposed,
    secondary structure, pLDDT confidence, proximity to functional sites)."""
    try:
        import structure_features as sf
    except Exception:
        return {"available": False}
    regions = [(int(r["start"]), int(r["end"])) for r in aging_regions]
    try:
        return sf.get_structural_context(uid, regions=regions)
    except Exception:
        return {"available": False}


def get_protein_context(uid: str) -> dict[str, Any]:
    """Assemble a structured, citable context for a protein.

    Returns a nested dict with sections: identity, function, domains, lipms,
    aging_signal, ml, orthologs, and sources (id -> url).
    """
    uid = uid.upper().strip()
    prot = _p("protein_features.parquet")
    row = prot[prot["uniprot_id"] == uid].drop_duplicates("uniprot_id")
    if len(row) == 0:
        return {"error": f"UniProt ID {uid} not found in the LiP-MS dataset."}
    r = row.iloc[0]

    peptides = _p("lipms_peptide_unified.parquet")
    uni = _uniprot_enrichment(uid)

    ctx: dict[str, Any] = {
        "identity": {
            "uniprot_id": uid,
            "gene_symbol": _safe(r.get("gene_symbol")),
            "recommended_name": _safe(uni.get("recommended_name")),
            "sequence_length": int(r["seq_length"]) if pd.notna(r.get("seq_length")) else None,
            "molecular_weight_kda": round(r["mol_weight"] / 1000, 1) if pd.notna(r.get("mol_weight")) else None,
        },
        "function": {
            "text": _safe(uni.get("function_text")),
            "go_biological_process": uni.get("go_bp", []),
            "go_molecular_function": uni.get("go_mf", []),
            "go_cellular_component": uni.get("go_cc", []),
            "keywords": uni.get("keywords", []),
            "subcellular_location": uni.get("subcellular_location", []),
        },
        "protein_annotations": {
            "log2_empai": round(float(r["log2_empai"]), 2) if pd.notna(r.get("log2_empai")) else None,
            "degree_centrality": round(float(r["degree_centrality"]), 4) if pd.notna(r.get("degree_centrality")) else None,
            "hsp90_client": _safe(r.get("hsp90_client"), None),
            "tric_client": _safe(r.get("tric_client"), None),
            "tissue_specificity": _safe(r.get("tissue_specificity")),
            "longevity": _safe(r.get("longevity"), None),
            "age_ribosome_pausing": _safe(r.get("age_ribosome_pausing"), None),
            "age_protein_turnover": _safe(r.get("age_protein_turnover"), None),
            "ubiquitination": _safe(r.get("ubiquitination"), None),
        },
        "domains": _domains(uid),
        "lipms": _peptide_summary(uid, peptides),
        "aging_signal": _lipms_aging_signal(uid, peptides),
    }
    # Annotate aging regions with the domains they overlap (using ALL InterPro entries, not just top-10 by size)
    ipr = _p("interpro_domains.parquet")
    all_domains = ipr[ipr["uniprot_id"] == uid]
    all_domains_list = [
        {"interpro_id": dr["interpro_id"], "name": dr["interpro_name"], "type": dr["domain_type"],
         "start": int(dr["domain_start"]), "end": int(dr["domain_end"])}
        for _, dr in all_domains.iterrows()
    ]
    ctx["aging_signal"]["top_regions"] = _annotate_regions_with_domains(
        ctx["aging_signal"]["top_regions"], all_domains_list)
    gene_for_url = r.get("gene_symbol") or uid
    ctx["variant_overlap"] = _variant_overlap(uid, ctx["aging_signal"]["top_regions"])
    ctx["structural_context"] = _structural_context(uid, ctx["aging_signal"]["top_regions"])
    ctx.update({
        "ml": _ml_predictions(uid),
        "orthologs": _orthologs(uid),
        "sources": {
            "S1": f"https://www.uniprot.org/uniprotkb/{uid}/entry",
            "S2": f"https://www.ebi.ac.uk/interpro/protein/UniProt/{uid}/",
            "S3": f"https://alphafold.ebi.ac.uk/entry/{uid}",
            "S4": f"https://wormbase.org/species/c_elegans/gene/{gene_for_url}",
            "S5": "This LiP-MS aging dataset (Sui et al.) — Table S1/S3, peptide-level conformational log2FC across 9 conditions.",
            "S6": "This ML predictor — LightGBM + ESM-2 35M PCA, isotonic-calibrated. See app/MODEL_CARD.md.",
            "S7": "ClinVar variants via UniProt Variation API, mapped worm→human via pairwise BLOSUM62 alignment. See VARIANTS_README.md.",
            "S8": "AlphaFold-derived structural features (biotite SASA/RSA, three-state secondary structure, pLDDT confidence) + CA–CA distance to UniProt functional sites. See app/STRUCTURE_README.md.",
        },
    })
    return ctx


if __name__ == "__main__":
    import sys
    uid = sys.argv[1] if len(sys.argv) > 1 else "P18948"
    ctx = get_protein_context(uid)
    print(json.dumps(ctx, indent=2, default=str))
