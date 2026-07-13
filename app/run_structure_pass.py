"""Batch per-residue structural features across all measured proteins with a local AF PDB."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import structure_features as sf

APP = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(APP, "data")

af = pd.read_parquet(os.path.join(DATA, "alphafold_index.parquet"))
per_res = pd.read_parquet(os.path.join(DATA, "lipms_per_residue_agg.parquet"))
measured = set(per_res["uniprot_id"].unique())

targets = af[af["uniprot_id"].isin(measured)].reset_index(drop=True)
print(f"targets: {len(targets)} measured proteins with local AF", flush=True)

def one(row):
    uid, path = row.uniprot_id, row.pdb_path
    try:
        df = sf.compute_structure_features(path)
        if not len(df):
            return uid, None, "empty"
        df.insert(0, "uniprot_id", uid)
        return uid, df, None
    except Exception as e:
        return uid, None, f"{type(e).__name__}: {e}"

results, errors = [], {}
t0 = time.time()
done = 0
with ThreadPoolExecutor(max_workers=6) as ex:
    futs = [ex.submit(one, row) for row in targets.itertuples(index=False)]
    for fut in as_completed(futs):
        uid, df, err = fut.result()
        done += 1
        if df is not None:
            results.append(df)
        else:
            errors[uid] = err
        if done % 500 == 0:
            print(f"  {done}/{len(targets)} @ {time.time()-t0:.0f}s  (ok={len(results)} err={len(errors)})", flush=True)

allf = pd.concat(results, ignore_index=True)
# compact dtypes
allf["residue"] = allf["residue"].astype("int32")
allf["sasa"] = allf["sasa"].astype("float32")
allf["rsa"] = allf["rsa"].astype("float32")
allf["plddt"] = allf["plddt"].astype("float32")
allf["sse"] = allf["sse"].astype("category")
out = os.path.join(DATA, "structural_features.parquet")
allf.to_parquet(out, index=False, compression="zstd")

print(f"\nDONE {time.time()-t0:.0f}s", flush=True)
print(f"proteins_ok={allf['uniprot_id'].nunique()} residues={len(allf)} errors={len(errors)}", flush=True)
print(f"pct_buried={100*allf['buried'].mean():.1f} pct_disordered={100*allf['disordered'].mean():.1f}", flush=True)
print(f"sse_dist={allf['sse'].value_counts(normalize=True).round(3).to_dict()}", flush=True)
print(f"plddt_median={allf['plddt'].median():.1f}", flush=True)
if errors:
    import json; json.dump(errors, open(os.path.join(DATA,"structure_pass_errors.json"),"w"))
    print(f"sample_errors={list(errors.items())[:5]}", flush=True)
print(f"wrote {out}", flush=True)
