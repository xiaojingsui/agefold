"""Per-residue 3D distance (CA-CA) to nearest functional site, from local AF structures."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import biotite.structure.io.pdb as pdbio
import biotite.structure as struc

APP = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(APP, "data")
NEAR_SITE_A = 8.0  # residue within 8 Å of a site CA → "near site"

af = pd.read_parquet(os.path.join(DATA, "alphafold_index.parquet")).set_index("uniprot_id")
sites = pd.read_parquet(os.path.join(DATA, "uniprot_sites.parquet"))
# Expand multi-residue site ranges into individual site residue numbers, keep type priority
TYPE_PRIORITY = {"Active site": 0, "Binding site": 1, "Site": 2, "Disulfide bond": 3}
site_res = {}  # uid -> list of (resnum, type)
for r in sites.itertuples(index=False):
    site_res.setdefault(r.uniprot_id, [])
    for pos in range(int(r.start), int(r.end) + 1):
        site_res[r.uniprot_id].append((pos, r.type))

def ca_coords(path):
    arr = pdbio.PDBFile.read(path).get_structure(model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    ca = arr[arr.atom_name == "CA"]
    return {int(rid): ca.coord[i] for i, rid in enumerate(ca.res_id)}

def one(uid):
    try:
        srs = site_res.get(uid)
        if not srs:
            return uid, None, "no_sites"
        coords = ca_coords(af.loc[uid, "pdb_path"])
        site_pts = [(rn, tp, coords[rn]) for rn, tp in srs if rn in coords]
        if not site_pts:
            return uid, None, "site_res_not_in_struct"
        pts = np.array([p[2] for p in site_pts])
        rows = []
        for rn, xyz in coords.items():
            d = np.linalg.norm(pts - xyz, axis=1)
            j = int(np.argmin(d))
            rows.append((uid, rn, site_pts[j][1], int(site_pts[j][0]), float(d[j]), bool(d[j] <= NEAR_SITE_A)))
        return uid, rows, None
    except Exception as e:
        return uid, None, f"{type(e).__name__}: {e}"

if __name__ == "__main__":
    uids = [u for u in site_res if u in af.index]
    print(f"computing site proximity for {len(uids)} proteins with sites", flush=True)
    all_rows, errs = [], {}
    t0 = time.time(); done = 0
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(one, u) for u in uids]
        for fut in as_completed(futs):
            uid, rows, err = fut.result()
            done += 1
            if rows: all_rows.extend(rows)
            elif err not in (None,): errs[uid] = err
            if done % 300 == 0:
                print(f"  {done}/{len(uids)} @ {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(all_rows, columns=["uniprot_id","residue","nearest_site_type","nearest_site_resnum","distance_ca","near_site"])
    df["residue"] = df["residue"].astype("int32")
    df["nearest_site_resnum"] = df["nearest_site_resnum"].astype("int32")
    df["distance_ca"] = df["distance_ca"].astype("float32")
    out = os.path.join(DATA, "site_proximity.parquet")
    df.to_parquet(out, index=False, compression="zstd")
    print(f"\nDONE {time.time()-t0:.0f}s proteins={df['uniprot_id'].nunique()} residues={len(df)} near_site_rows={int(df['near_site'].sum())} errors={len(errs)}", flush=True)
    print(f"wrote {out}", flush=True)
