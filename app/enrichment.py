"""
UniProt + InterPro live lookups with SQLite-backed on-disk cache.
Falls back to local InterPro parquet when offline.
"""
from __future__ import annotations
import json, sqlite3, time
from pathlib import Path
from typing import Optional

import requests
import pandas as pd

CACHE_PATH = Path(__file__).parent / "data" / "enrichment_cache.sqlite"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{uid}.json"
INTERPRO_URL = "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{uid}/"


def _ensure_cache() -> sqlite3.Connection:
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(CACHE_PATH))
    conn.execute("""
        CREATE TABLE IF NOT EXISTS enrichment (
            uniprot_id TEXT NOT NULL,
            source     TEXT NOT NULL,
            payload    TEXT NOT NULL,
            fetched_at REAL NOT NULL,
            PRIMARY KEY (uniprot_id, source)
        )
    """)
    conn.commit()
    return conn


def _cache_get(conn, uid: str, source: str) -> Optional[dict]:
    row = conn.execute("SELECT payload FROM enrichment WHERE uniprot_id=? AND source=?",
                       (uid, source)).fetchone()
    return json.loads(row[0]) if row else None


def _cache_put(conn, uid: str, source: str, payload: dict) -> None:
    conn.execute("INSERT OR REPLACE INTO enrichment(uniprot_id, source, payload, fetched_at) VALUES (?,?,?,?)",
                 (uid, source, json.dumps(payload), time.time()))
    conn.commit()


# ---------------- UniProt ----------------
def _parse_uniprot(raw: dict) -> dict:
    """Extract the fields the protein card needs from a UniProt JSON blob."""
    out = {"function": None, "recommended_name": None,
           "go_bp": [], "go_mf": [], "go_cc": [],
           "keywords": [], "subcellular_location": []}
    # recommended name
    rn = raw.get("proteinDescription", {}).get("recommendedName", {})
    out["recommended_name"] = rn.get("fullName", {}).get("value")
    # function
    for c in raw.get("comments", []):
        ct = c.get("commentType")
        if ct == "FUNCTION":
            for t in c.get("texts", []):
                v = t.get("value")
                if v:
                    out["function"] = (out["function"] + " " + v) if out["function"] else v
        elif ct == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations", []):
                v = loc.get("location", {}).get("value")
                if v: out["subcellular_location"].append(v)
    # keywords
    for kw in raw.get("keywords", []):
        v = kw.get("name")
        if v: out["keywords"].append(v)
    # GO by namespace
    for x in raw.get("uniProtKBCrossReferences", []):
        if x.get("database") != "GO":
            continue
        term_ns = None; term_name = None
        for p in x.get("properties", []):
            if p.get("key") == "GoTerm":
                v = p.get("value", "")
                if ":" in v:
                    term_ns, term_name = v.split(":", 1)
        entry = {"id": x.get("id"), "name": term_name}
        if term_ns == "P": out["go_bp"].append(entry)
        elif term_ns == "F": out["go_mf"].append(entry)
        elif term_ns == "C": out["go_cc"].append(entry)
    return out


def fetch_uniprot(uid: str, force: bool = False) -> dict:
    """Return parsed UniProt enrichment for a protein; cached."""
    conn = _ensure_cache()
    if not force:
        hit = _cache_get(conn, uid, "uniprot")
        if hit is not None:
            return hit
    try:
        r = requests.get(UNIPROT_URL.format(uid=uid), timeout=15,
                         headers={"Accept": "application/json"})
        if r.status_code != 200:
            payload = {"error": f"UniProt {r.status_code}", "recommended_name": None,
                       "function": None, "go_bp": [], "go_mf": [], "go_cc": [],
                       "keywords": [], "subcellular_location": []}
        else:
            payload = _parse_uniprot(r.json())
    except Exception as e:
        payload = {"error": str(e), "recommended_name": None, "function": None,
                   "go_bp": [], "go_mf": [], "go_cc": [], "keywords": [],
                   "subcellular_location": []}
    _cache_put(conn, uid, "uniprot", payload)
    return payload


# ---------------- InterPro (live) ----------------
def _parse_interpro(raw: dict) -> dict:
    entries = []
    for res in raw.get("results", []):
        meta = res.get("metadata", {})
        # locations may live in res['proteins'][*]['entry_protein_locations']
        locs = []
        for prot in res.get("proteins", []) or []:
            for loc in prot.get("entry_protein_locations", []) or []:
                for frag in loc.get("fragments", []) or []:
                    locs.append((int(frag.get("start", 0)), int(frag.get("end", 0))))
        entries.append({
            "accession": meta.get("accession"),
            "name": meta.get("name"),
            "type": meta.get("type"),
            "source_database": meta.get("source_database"),
            "locations": locs,
        })
    return {"count": raw.get("count", len(entries)), "entries": entries}


def fetch_interpro_live(uid: str, force: bool = False) -> dict:
    conn = _ensure_cache()
    if not force:
        hit = _cache_get(conn, uid, "interpro")
        if hit is not None:
            return hit
    try:
        r = requests.get(INTERPRO_URL.format(uid=uid), timeout=15,
                         headers={"Accept": "application/json"})
        payload = _parse_interpro(r.json()) if r.status_code == 200 else {"error": f"InterPro {r.status_code}", "count": 0, "entries": []}
    except Exception as e:
        payload = {"error": str(e), "count": 0, "entries": []}
    _cache_put(conn, uid, "interpro", payload)
    return payload


# ---------------- HGNC -> Human UniProt ----------------
def fetch_hgnc_to_uniprot(hgnc_id: str, force: bool = False) -> dict:
    """
    Given an HGNC ID like 'HGNC:5241', return the reviewed human UniProt entry:
    {accession, gene_name, protein_name, length, error?}.
    """
    key = hgnc_id.replace(":", "_")
    conn = _ensure_cache()
    if not force:
        hit = _cache_get(conn, key, "hgnc")
        if hit is not None:
            return hit
    hgnc_num = hgnc_id.split(":")[-1]
    try:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/search",
            params={
                "query": f"xref:hgnc-{hgnc_num} AND organism_id:9606 AND reviewed:true",
                "fields": "accession,id,gene_names,protein_name,length",
                "format": "json", "size": 1,
            },
            timeout=15,
        )
        if r.status_code != 200 or not r.json().get("results"):
            payload = {"accession": None, "gene_name": None, "protein_name": None, "length": None,
                       "error": f"UniProt {r.status_code}" if r.status_code != 200 else "no hits"}
        else:
            hit = r.json()["results"][0]
            gene = hit.get("genes", [{}])
            gene_name = gene[0].get("geneName", {}).get("value") if gene else None
            payload = {
                "accession": hit.get("primaryAccession"),
                "gene_name": gene_name,
                "protein_name": hit.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                "length": hit.get("sequence", {}).get("length"),
            }
    except Exception as e:
        payload = {"accession": None, "gene_name": None, "protein_name": None, "length": None, "error": str(e)}
    _cache_put(conn, key, "hgnc", payload)
    return payload


# ---------------- InterPro fallback (local parquet) ----------------
def interpro_local(uid: str, local_df: pd.DataFrame) -> list[dict]:
    sub = local_df[local_df["uniprot_id"] == uid]
    return [
        {"accession": r["interpro_id"], "name": r["interpro_name"],
         "type": r["domain_type"], "source_database": "local",
         "locations": [(int(r["domain_start"]), int(r["domain_end"]))]}
        for _, r in sub.iterrows()
    ]
