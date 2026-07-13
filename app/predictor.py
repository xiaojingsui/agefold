"""
Inference helper for the residue-level aging-destabilization classifier.

Given a UniProt ID, load its features from the parquet layer and return a
per-residue predicted probability (isotonic-calibrated).
"""
from __future__ import annotations
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
import lightgbm as lgb

APP_DIR = Path(__file__).parent
DATA_DIR = APP_DIR / "data"

# Feature order must match training
FEATURES = [
    "is_semi","is_half","is_full","pk_hydrophobic",
    "avg_residue_depth","plddt","is_idr",
    "ss_helix","ss_sheet","ss_turn","ss_coil",
    "has_domain","in_active_site","in_binding_site",
    "rel_pos","n_peptides",
    "seq_length","mol_weight","charge","log2_empai","degree_centrality",
    "has_hsp90","has_tric","tissue_specific",
    "is_prolongevity","is_antilongevity","ribosome_pausing","is_ubiquitinated",
]

_state: dict = {}


def _load():
    if "model" in _state:
        return _state
    _state["model"] = lgb.Booster(model_file=str(DATA_DIR / "model_full.txt"))
    with open(DATA_DIR / "isotonic_calibrator.pkl", "rb") as f:
        _state["iso"] = pickle.load(f)
    _state["preds"] = pd.read_parquet(DATA_DIR / "residue_predictions.parquet")
    return _state


def score_uniprot(uniprot_id: str) -> pd.DataFrame:
    """
    Return DataFrame(residue, p_destabilized, observed_in_lipms) for a protein.
    Uses the precomputed table for observed residues (fast). For proteins not
    in the training set, returns an empty DataFrame — use `score_from_matrix`
    with a fresh feature matrix.
    """
    s = _load()
    return s["preds"].query("uniprot_id == @uniprot_id").sort_values("residue").reset_index(drop=True)


def score_from_matrix(X: pd.DataFrame) -> np.ndarray:
    """Apply model + isotonic calibrator to a feature matrix (columns = FEATURES)."""
    s = _load()
    p_raw = s["model"].predict(X[FEATURES].values)
    return s["iso"].predict(p_raw).astype("float32")


def top_residues(uniprot_id: str, k: int = 10) -> pd.DataFrame:
    """Return the top-k residues by predicted probability for a protein."""
    d = score_uniprot(uniprot_id)
    return d.nlargest(k, "p_destabilized")
