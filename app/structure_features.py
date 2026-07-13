"""
structure_features.py — per-residue structural context from local AlphaFold PDBs.

Computes, for each residue of a protein:
  - sasa      : absolute solvent-accessible surface area (Å², Shrake-Rupley)
  - rsa       : relative solvent accessibility = sasa / max-ASA(residue type)
  - buried    : rsa < BURIED_RSA  (True = buried core, False = exposed)
  - sse       : three-state secondary structure  'H' helix / 'E' sheet / 'C' coil
  - plddt     : AlphaFold per-residue confidence (mean CA B-factor)
  - disordered: plddt < DISORDER_PLDDT (very low confidence → likely disordered)

pLDDT note: AlphaFold PDBs store per-residue pLDDT in the B-factor column
(0-100). We take the CA atom's value per residue.

Residue numbering is 1-based UniProt numbering, which matches the local
AlphaFold PDBs 1:1 (established for this dataset — no offset).
"""
from __future__ import annotations
import numpy as np
import pandas as pd

import biotite.structure as struc
import biotite.structure.io.pdb as pdbio

# --- thresholds ---
BURIED_RSA = 0.20        # < 20% relative accessibility → buried core
DISORDER_PLDDT = 50.0    # AlphaFold pLDDT < 50 → very low confidence / likely disordered

# Tien et al. 2013 theoretical maximum accessible surface area (Å²) per residue.
# (Used to normalize SASA into relative accessibility.)
MAX_ASA_TIEN = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
_SSE_MAP = {"a": "H", "b": "E", "c": "C"}  # biotite annotate_sse: a=alpha, b=beta, c=coil


def compute_structure_features(pdb_path: str) -> pd.DataFrame:
    """Per-residue structural features for one AlphaFold PDB.
    Returns DataFrame indexed by residue number (1-based)."""
    pdb_file = pdbio.PDBFile.read(pdb_path)
    # B-factor carries pLDDT for AlphaFold models
    arr = pdb_file.get_structure(model=1, extra_fields=["b_factor"])
    # protein atoms only, first (only) chain
    arr = arr[struc.filter_amino_acids(arr)]
    if arr.array_length() == 0:
        return pd.DataFrame(columns=["residue", "res_name", "sasa", "rsa", "buried", "sse", "plddt", "disordered"])

    # --- SASA (per atom) → sum per residue ---
    atom_sasa = struc.sasa(arr, vdw_radii="Single")
    res_ids = arr.res_id
    res_names = arr.res_name

    # aggregate per residue
    uniq_res = np.unique(res_ids)
    sasa_per_res, name_per_res, plddt_per_res = [], [], []
    for r in uniq_res:
        mask = res_ids == r
        s = np.nansum(atom_sasa[mask])
        sasa_per_res.append(s)
        name_per_res.append(res_names[mask][0])
        # pLDDT = CA b_factor if present, else mean
        ca = mask & (arr.atom_name == "CA")
        plddt_per_res.append(float(arr.b_factor[ca][0]) if ca.any() else float(np.nanmean(arr.b_factor[mask])))

    df = pd.DataFrame({
        "residue": uniq_res.astype(int),
        "res_name": name_per_res,
        "sasa": np.array(sasa_per_res, dtype=float),
        "plddt": np.array(plddt_per_res, dtype=float),
    })
    df["rsa"] = [min(s / MAX_ASA_TIEN.get(n, np.nan), 1.5) if MAX_ASA_TIEN.get(n) else np.nan
                 for s, n in zip(df["sasa"], df["res_name"])]
    df["buried"] = df["rsa"] < BURIED_RSA
    df["disordered"] = df["plddt"] < DISORDER_PLDDT

    # --- secondary structure (three-state) ---
    try:
        sse = struc.annotate_sse(arr)  # one letter per residue, in residue order
        # annotate_sse returns per-residue array aligned to unique residues in structure order
        sse3 = np.array([_SSE_MAP.get(x, "C") for x in sse])
        if len(sse3) == len(df):
            df["sse"] = sse3
        else:
            # length mismatch (rare) → align by index up to min length
            out = np.full(len(df), "C", dtype=object)
            out[:min(len(df), len(sse3))] = sse3[:min(len(df), len(sse3))]
            df["sse"] = out
    except Exception:
        df["sse"] = "C"

    return df[["residue", "res_name", "sasa", "rsa", "buried", "sse", "plddt", "disordered"]]


# ------------------------------------------------------------------
# Precomputed lookups (loaded lazily; used by the app + chat)
# ------------------------------------------------------------------
import os
from pathlib import Path

_DATA = Path(__file__).parent / "data"
_cache = {}

def _load(name):
    if name not in _cache:
        p = _DATA / name
        _cache[name] = pd.read_parquet(p) if p.exists() else None
    return _cache[name]


_SSE_NAME = {"H": "helix", "E": "sheet/strand", "C": "loop/coil"}


def _interpret(buried_frac, dom_sse, mean_plddt, near_site, site_type, disorder_frac):
    """One-line plain structural interpretation of an aging region."""
    # confidence gate first
    if mean_plddt < DISORDER_PLDDT or disorder_frac > 0.5:
        return ("This region is low-confidence in the AlphaFold model (likely flexible or "
                "disordered) — shape changes here are harder to interpret structurally.")
    core = buried_frac >= 0.5
    loc = "buried in the folded core" if core else "on the solvent-exposed surface"
    sse_txt = _SSE_NAME.get(dom_sse, "loop/coil")
    parts = [f"Mostly {loc}, predominantly {sse_txt}"]
    if near_site:
        st = {"Binding site": "a ligand/metal-binding site", "Active site": "the active site",
              "Disulfide bond": "a disulfide bond", "Site": "an annotated functional site"}.get(site_type, "a functional site")
        parts.append(f"and lies close to {st}")
        if core:
            tail = "Destabilization of a buried, functionally coupled region like this is the kind most likely to affect fold stability or activity."
        else:
            tail = "A change next to a functional site, even at the surface, can affect binding or regulation."
    else:
        if core:
            tail = "Buried-core changes are more likely to affect overall fold stability than surface changes."
        else:
            tail = "Surface changes are often better tolerated and may reflect altered interactions rather than loss of fold."
    return ". ".join([", ".join(parts), tail])


def get_structural_context(uid, regions=None):
    """Structural context for a protein, optionally summarized over aging regions.

    regions: list of (start, end) 1-based inclusive. If None, only per-protein
    availability is returned.
    Returns dict: {available, per_protein:{...}, regions:[{start,end,pct_buried,
      dominant_sse, mean_plddt, low_confidence, near_site, nearest_site_type,
      min_site_distance, interpretation}]}
    """
    sf = _load("structural_features.parquet")
    if sf is None:
        return {"available": False}
    s = sf[sf["uniprot_id"] == uid]
    if not len(s):
        return {"available": False}
    prox = _load("site_proximity.parquet")
    p = prox[prox["uniprot_id"] == uid] if prox is not None else None
    p_idx = p.set_index("residue") if p is not None and len(p) else None

    out = {"available": True,
           "per_protein": {
               "n_residues": int(len(s)),
               "pct_buried": round(100 * s["buried"].mean(), 1),
               "median_plddt": round(float(s["plddt"].median()), 1),
               "has_sites": bool(p_idx is not None),
           },
           "regions": []}
    if not regions:
        return out

    s_idx = s.set_index("residue")
    for (start, end) in regions:
        seg = s_idx[(s_idx.index >= start) & (s_idx.index <= end)]
        if not len(seg):
            continue
        buried_frac = float(seg["buried"].mean())
        disorder_frac = float(seg["disordered"].mean())
        dom_sse = seg["sse"].mode().iloc[0] if len(seg["sse"].mode()) else "C"
        mean_plddt = float(seg["plddt"].mean())
        near_site = False; nearest_type = None; min_dist = None
        if p_idx is not None:
            pseg = p_idx[(p_idx.index >= start) & (p_idx.index <= end)]
            if len(pseg):
                near_site = bool(pseg["near_site"].any())
                min_dist = round(float(pseg["distance_ca"].min()), 1)
                nearest_type = pseg.loc[pseg["distance_ca"].idxmin(), "nearest_site_type"]
        out["regions"].append({
            "start": int(start), "end": int(end),
            "pct_buried": round(100 * buried_frac, 1),
            "dominant_sse": _SSE_NAME.get(dom_sse, "loop/coil"),
            "mean_plddt": round(mean_plddt, 1),
            "low_confidence": bool(mean_plddt < DISORDER_PLDDT or disorder_frac > 0.5),
            "near_site": near_site,
            "nearest_site_type": nearest_type,
            "min_site_distance": min_dist,
            "interpretation": _interpret(buried_frac, dom_sse, mean_plddt, near_site, nearest_type, disorder_frac),
        })
    return out
