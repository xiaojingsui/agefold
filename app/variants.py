"""
Variant retrieval + worm↔human coordinate mapping.

Provides:
- fetch_human_variants(human_uid): UniProt Variation API, cached in SQLite
- align_worm_human(worm_seq, human_seq): pairwise BLOSUM62 alignment
- map_variants_to_worm(variants, alignment, human_uid, worm_uid): apply mapping
"""
from __future__ import annotations
import json, sqlite3, time
from pathlib import Path
from typing import Any, Optional
import requests

APP_DIR = Path(__file__).parent
DATA_DIR = APP_DIR / "data"
CACHE_PATH = DATA_DIR / "enrichment_cache.sqlite"

VARIATION_URL = "https://www.ebi.ac.uk/proteins/api/variation/{uid}"
TIMEOUT = 20

# Clinical significance labels we consider "pathogenic" (germline disease-linked)
PATHOGENIC_LABELS = {
    "Likely pathogenic", "Pathogenic",
    "Likely pathogenic, low penetrance", "Pathogenic, low penetrance",
    "Disease", "Disease causing", "Association",
}

def _cache_get(source: str, key: str) -> Optional[dict]:
    if not CACHE_PATH.exists(): return None
    conn = sqlite3.connect(str(CACHE_PATH))
    try:
        row = conn.execute("SELECT payload FROM enrichment WHERE uniprot_id=? AND source=?",
                            (key, source)).fetchone()
        return json.loads(row[0]) if row else None
    finally:
        conn.close()

def _cache_put(source: str, key: str, payload: dict) -> None:
    # Self-create the DB + table so a variants-first write on a fresh deploy
    # (where enrichment_cache.sqlite is gitignored and absent) does not raise
    # "no such table: enrichment". Schema matches enrichment.py's table.
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(CACHE_PATH))
    try:
        conn.execute("""
            CREATE TABLE IF NOT EXISTS enrichment (
                uniprot_id TEXT NOT NULL,
                source     TEXT NOT NULL,
                payload    TEXT NOT NULL,
                fetched_at REAL NOT NULL,
                PRIMARY KEY (uniprot_id, source)
            )
        """)
        conn.execute(
            "INSERT OR REPLACE INTO enrichment (uniprot_id, source, payload, fetched_at) VALUES (?,?,?,?)",
            (key, source, json.dumps(payload), time.time()))
        conn.commit()
    finally:
        conn.close()


def _parse_variation_response(raw: dict) -> dict:
    """Extract features → normalized variant dicts."""
    features = raw.get("features") or []
    variants = []
    for f in features:
        # Skip non-variant features
        if f.get("type") not in ("VARIANT", "Variant"):
            continue
        # Position
        try:
            begin = int(f.get("begin", 0))
            end = int(f.get("end", begin))
        except (ValueError, TypeError):
            continue
        if begin <= 0: continue

        # Collect clinical significance labels
        cs = f.get("clinicalSignificances") or []
        cs_types = []
        for entry in cs:
            t = entry.get("type")
            if t: cs_types.append(t)

        # Associations → diseases
        assoc = f.get("association") or []
        diseases = []
        for a in assoc:
            d = a.get("name") or a.get("description")
            if d: diseases.append(d)

        variants.append({
            "position": begin,
            "end": end,
            "wildType": f.get("wildType"),
            "mutatedType": f.get("mutatedType") or f.get("alternativeSequence"),
            "clinicalSignificances": cs_types,
            "is_pathogenic": any(t in PATHOGENIC_LABELS for t in cs_types),
            "somaticStatus": f.get("somaticStatus"),
            "diseases": diseases[:5],  # cap
            "consequence": f.get("consequenceType"),
            "sources": f.get("sourceType") or "",
        })
    return {
        "accession": raw.get("accession"),
        "sequence": raw.get("sequence"),
        "sequence_length": len(raw.get("sequence") or ""),
        "n_variants": len(variants),
        "n_pathogenic": sum(1 for v in variants if v["is_pathogenic"]),
        "variants": variants,
    }


def fetch_human_variants(human_uid: str, force: bool = False,
                          session: Optional[requests.Session] = None) -> dict:
    """Return parsed variant dict for a human UniProt accession, using SQLite cache."""
    if not force:
        cached = _cache_get("uniprot_variation", human_uid)
        if cached is not None:
            return cached
    sess = session or requests
    url = VARIATION_URL.format(uid=human_uid)
    r = sess.get(url, headers={"Accept": "application/json"}, timeout=TIMEOUT)
    if r.status_code == 404:
        payload = {"accession": human_uid, "sequence": None, "sequence_length": 0,
                    "n_variants": 0, "n_pathogenic": 0, "variants": [], "error": "not_found"}
    elif r.status_code >= 400:
        payload = {"accession": human_uid, "n_variants": 0, "n_pathogenic": 0,
                    "variants": [], "error": f"http_{r.status_code}"}
    else:
        payload = _parse_variation_response(r.json())
    _cache_put("uniprot_variation", human_uid, payload)
    return payload


def align_worm_human(worm_seq: str, human_seq: str) -> dict:
    """Pairwise BLOSUM62 alignment. Returns {alignment_str, identity, human_to_worm mapping}."""
    from Bio import Align
    from Bio.Align import substitution_matrices
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"
    alns = aligner.align(worm_seq, human_seq)
    aln = alns[0]
    # Build human_pos → worm_pos (1-based); use the aligned coordinate arrays
    # aln.aligned is ((worm_starts, worm_ends), (human_starts, human_ends)) as arrays of blocks
    worm_blocks, human_blocks = aln.aligned
    human_to_worm = {}
    n_match = 0
    for (ws, we), (hs, he) in zip(worm_blocks, human_blocks):
        for offset in range(we - ws):
            worm_i = ws + offset  # 0-based
            human_i = hs + offset
            if worm_i < len(worm_seq) and human_i < len(human_seq):
                # 1-based mapping
                human_to_worm[human_i + 1] = worm_i + 1
                if worm_seq[worm_i] == human_seq[human_i]:
                    n_match += 1
    aln_len = sum(we - ws for (ws, we), _ in zip(worm_blocks, human_blocks))
    identity = n_match / aln_len if aln_len else 0.0
    return {
        "human_to_worm": human_to_worm,
        "n_aligned": aln_len,
        "identity": identity,
        "n_match": n_match,
        "score": float(aln.score),
    }


def map_variants_to_worm(variant_payload: dict, alignment: dict,
                          human_uid: str, worm_uid: str) -> list[dict]:
    """Apply the alignment to each variant, dropping ones that don't map."""
    mapped = []
    h2w = alignment["human_to_worm"]
    id_ = alignment["identity"]
    for v in variant_payload["variants"]:
        wp = h2w.get(v["position"])
        if wp is None: continue
        mapped.append({
            "worm_uniprot": worm_uid,
            "worm_residue": int(wp),
            "human_uid": human_uid,
            "human_pos": int(v["position"]),
            "wildtype_hs": v["wildType"],
            "mutant_hs": v["mutatedType"],
            "is_pathogenic": v["is_pathogenic"],
            "clinical": ";".join(v["clinicalSignificances"]),
            "diseases": v["diseases"],
            "alignment_identity": id_,
        })
    return mapped


if __name__ == "__main__":
    import sys
    uid = sys.argv[1] if len(sys.argv) > 1 else "P11142"
    v = fetch_human_variants(uid)
    print(json.dumps({"accession": v.get("accession"),
                       "n_variants": v.get("n_variants"),
                       "n_pathogenic": v.get("n_pathogenic"),
                       "sample": v.get("variants", [])[:3]}, indent=2))
