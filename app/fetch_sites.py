"""Fetch UniProt functional-site features (active/binding/metal/disulfide) for measured proteins.

Caches into enrichment_cache.sqlite under source='uniprot_sites' as JSON list of
{type, start, end, description, ligand}. Threaded with backoff.
"""
import os, sys, json, time, sqlite3, threading
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

APP = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(APP, "data")
CACHE = os.path.join(DATA, "enrichment_cache.sqlite")
URL = "https://rest.uniprot.org/uniprotkb/{uid}.json?fields=ft_act_site,ft_binding,ft_site,ft_disulfid"

SITE_TYPES = {"Active site", "Binding site", "Site", "Metal binding", "Disulfide bond"}
_local = threading.local()

def _conn():
    if not hasattr(_local, "c"):
        # Self-create the DB + table so a sites read/write on a fresh deploy
        # (where enrichment_cache.sqlite is gitignored and absent) does not
        # raise "no such table: enrichment". Schema matches enrichment.py.
        os.makedirs(DATA, exist_ok=True)
        _local.c = sqlite3.connect(CACHE, check_same_thread=False, timeout=30)
        _local.c.execute("""
            CREATE TABLE IF NOT EXISTS enrichment (
                uniprot_id TEXT NOT NULL,
                source     TEXT NOT NULL,
                payload    TEXT NOT NULL,
                fetched_at REAL NOT NULL,
                PRIMARY KEY (uniprot_id, source)
            )
        """)
        _local.c.commit()
    return _local.c

def _cached(uid):
    r = _conn().execute("SELECT payload FROM enrichment WHERE source='uniprot_sites' AND uniprot_id=?", (uid,)).fetchone()
    return json.loads(r[0]) if r else None

def _store(uid, sites):
    c = _conn()
    c.execute("INSERT OR REPLACE INTO enrichment(uniprot_id, source, payload, fetched_at) VALUES (?,?,?,?)",
              (uid, "uniprot_sites", json.dumps(sites), int(time.time())))
    c.commit()

def _parse(d):
    out = []
    for f in d.get("features", []):
        if f["type"] not in SITE_TYPES:
            continue
        loc = f.get("location", {})
        s = loc.get("start", {}).get("value"); e = loc.get("end", {}).get("value")
        if s is None:
            continue
        out.append({"type": f["type"], "start": int(s), "end": int(e) if e else int(s),
                    "description": f.get("description", ""),
                    "ligand": f.get("ligand", {}).get("name", "")})
    return out

def fetch_one(uid, retries=3):
    c = _cached(uid)
    if c is not None:
        return uid, c, "cached"
    for attempt in range(retries):
        try:
            r = requests.get(URL.format(uid=uid), timeout=20)
            if r.status_code == 200:
                sites = _parse(r.json())
                _store(uid, sites)
                return uid, sites, "fetched"
            if r.status_code == 429:
                time.sleep(2 ** attempt); continue
            return uid, [], f"http{r.status_code}"
        except Exception as e:
            if attempt == retries - 1:
                return uid, None, f"{type(e).__name__}"
            time.sleep(2 ** attempt)
    return uid, None, "retries"

if __name__ == "__main__":
    per_res = pd.read_parquet(os.path.join(DATA, "lipms_per_residue_agg.parquet"))
    uids = sorted(per_res["uniprot_id"].unique())
    print(f"fetching sites for {len(uids)} measured proteins", flush=True)
    t0 = time.time(); done = 0; with_sites = 0; errs = 0
    rows = []
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(fetch_one, u) for u in uids]
        for fut in as_completed(futs):
            uid, sites, status = fut.result()
            done += 1
            if sites is None:
                errs += 1
            else:
                if sites:
                    with_sites += 1
                    for s in sites:
                        rows.append({"uniprot_id": uid, **s})
            if done % 1000 == 0:
                print(f"  {done}/{len(uids)} @ {time.time()-t0:.0f}s (with_sites={with_sites} err={errs})", flush=True)
    df = pd.DataFrame(rows)
    out = os.path.join(DATA, "uniprot_sites.parquet")
    df.to_parquet(out, index=False)
    print(f"\nDONE {time.time()-t0:.0f}s proteins_with_sites={with_sites} site_rows={len(df)} errors={errs}", flush=True)
    if len(df):
        print(f"site_type_dist={df['type'].value_counts().to_dict()}", flush=True)
    print(f"wrote {out}", flush=True)
