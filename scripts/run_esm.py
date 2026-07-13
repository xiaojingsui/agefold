"""
Extract ESM-2 35M per-residue embeddings for all proteins in the training set.
Writes one parquet per protein under data/esm_chunks/<uid>.parquet.
Sequences longer than 1022 are chunked with 128-aa overlap; overlapping positions
are averaged.
"""
import os, json, time, sys
import numpy as np
import pandas as pd
import torch
import esm

torch.set_num_threads(max(1, os.cpu_count()-1))

DEVICE = "mps" if torch.backends.mps.is_available() else "cpu"
print(f"Device: {DEVICE}, torch threads: {torch.get_num_threads()}", flush=True)

print("Loading esm2_t12_35M_UR50D...", flush=True)
t0 = time.time()
model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
model = model.eval()
try:
    model = model.to(DEVICE)
except Exception as e:
    print(f"WARN: could not move to {DEVICE}: {e}; falling back to cpu", flush=True)
    DEVICE = "cpu"
bc = alphabet.get_batch_converter()
N_LAYER = 12
D = 480
print(f"Model loaded in {time.time()-t0:.1f}s", flush=True)

with open('data/protein_sequences.json') as f:
    seqs = json.load(f)
print(f"{len(seqs):,} sequences to embed", flush=True)

CHUNK = 1022
OVERLAP = 128

def _chunks(seq):
    n = len(seq)
    if n <= CHUNK:
        yield 0, seq
        return
    step = CHUNK - OVERLAP
    for start in range(0, n, step):
        end = min(start + CHUNK, n)
        yield start, seq[start:end]
        if end == n: break

def embed_one(uid, sd):
    seq = sd['seq']; first_res = sd['first_res']
    n = len(seq)
    if n == 0: return None
    accum = np.zeros((n, D), dtype=np.float32)
    counts = np.zeros(n, dtype=np.int32)
    for start, chunk in _chunks(seq):
        try:
            _, _, toks = bc([(uid, chunk)])
            toks = toks.to(DEVICE)
            with torch.no_grad():
                out = model(toks, repr_layers=[N_LAYER])
            emb = out["representations"][N_LAYER][0, 1:-1].cpu().numpy().astype(np.float32)
        except Exception as e:
            print(f"FAIL {uid} chunk@{start}: {e}", flush=True)
            return None
        accum[start:start+len(chunk)] += emb
        counts[start:start+len(chunk)] += 1
    accum /= np.maximum(counts[:, None], 1)
    return pd.DataFrame({
        'uniprot_id': uid,
        'residue': np.arange(first_res, first_res + n),
        **{f'e{i}': accum[:, i] for i in range(D)}
    })

os.makedirs('data/esm_chunks', exist_ok=True)
uids = sorted(seqs.keys())
done_files = set(os.listdir('data/esm_chunks'))
todo = [u for u in uids if f'{u}.parquet' not in done_files]
print(f"To embed: {len(todo):,} (skipping {len(uids)-len(todo)} already done)", flush=True)

t0 = time.time()
for i, u in enumerate(todo):
    df = embed_one(u, seqs[u])
    if df is None: continue
    df.to_parquet(f'data/esm_chunks/{u}.parquet', index=False, compression='zstd')
    if (i+1) % 100 == 0 or i == len(todo)-1:
        el = time.time() - t0
        rate = (i+1) / max(el, 1)
        eta = (len(todo) - i - 1) / max(rate, 1e-6)
        print(f"  {i+1:,}/{len(todo):,}  ({el:.0f}s, {rate:.2f}/s, ETA {eta/60:.1f} min)", flush=True)

print(f"\nAll done in {time.time()-t0:.0f}s", flush=True)
