# Deploying AgeFold to Streamlit Community Cloud

A free, permanent public URL — the app runs on Streamlit's servers, no laptop
required. This is the recommended path for linking from a paper or sharing
long-term.

## What makes this deployable

Two things were converted so the app runs off your Mac:

1. **AlphaFold structures are fetched on demand.** The app previously read
   19,694 PDBs from local OneDrive paths that don't exist on a server. It now
   calls `app_core.get_pdb_text(uid)`, which uses the local file if present
   (fast path on your Mac) and otherwise queries the EBI AlphaFold API
   (`alphafold.ebi.ac.uk/api/prediction/<uid>`) for the *current* model URL,
   downloads it, and caches it under `app/data/af_cache/`. The bundled model
   version drifts (v2 → v6 over time), so the URL is always read from the API,
   never hardcoded.

2. **The 512 MB enrichment cache is not committed.**
   `app/data/enrichment_cache.sqlite` is gitignored. All three writers
   (`enrichment.py`, `variants.py`, `fetch_sites.py`) self-create the DB and
   the `enrichment` table on first use (`CREATE TABLE IF NOT EXISTS` + `mkdir`),
   and every read guards on the file/table being absent. So on a fresh server
   the cache rebuilds lazily from live UniProt / InterPro / ClinVar fetches —
   the first view of each protein is a little slower, then it's cached.

The committed repo is ~128 MB (parquets + code + assets); largest single file
is 36 MB — comfortably within GitHub and Community Cloud limits.

## What is NOT in the repo (by design)

`.gitignore` excludes:
- `.streamlit/secrets.toml` — your Anthropic API key. **Never commit this.**
- `app/data/enrichment_cache.sqlite` — the 512 MB cache (rebuilds lazily).
- `app/data/af_cache/` — on-demand-fetched structures (rebuild lazily).
- `__pycache__/`, `*.pyc`, `.DS_Store`, editor dirs.

---

## Step 1 — Initialize git and commit (in your Terminal)

Git can't be run from the assistant sandbox — do this in Terminal on your Mac:

```bash
cd ~/Documents/lipms-viewer
git init
git add .
git status      # CONFIRM: no secrets.toml, no enrichment_cache.sqlite, no af_cache/
git commit -m "AgeFold: conversational AI viewer for C. elegans LiP-MS aging data"
```

If `git status` lists `.streamlit/secrets.toml` or
`app/data/enrichment_cache.sqlite`, stop — the `.gitignore` isn't being picked
up. (It lives at the repo root next to this file.)

## Step 2 — Create the GitHub repo and push

**Option A — GitHub CLI** (if you have `gh`):
```bash
gh repo create agefold --public --source=. --remote=origin --push
```

**Option B — github.com:**
1. github.com → **New repository** → name `agefold` → **Public** →
   do **not** add a README/.gitignore/license → **Create repository**.
2. Back in Terminal:
   ```bash
   git remote add origin https://github.com/YOUR_USERNAME/agefold.git
   git branch -M main
   git push -u origin main
   ```

## Step 3 — Deploy on Streamlit Community Cloud

1. Go to **https://share.streamlit.io** and sign in with GitHub.
2. **Create app → Deploy a public app from GitHub**.
3. Fill in:
   - **Repository:** `YOUR_USERNAME/agefold`
   - **Branch:** `main`
   - **Main file path:** `app/streamlit_app.py`   ← important, the app is in `app/`
4. Open **Advanced settings → Secrets** and paste your key (TOML format):
   ```toml
   anthropic_api_key = "sk-ant-...your real key..."
   ```
   This is the same value as your local `.streamlit/secrets.toml`, which is not
   in the repo. Without it the app still runs — chat falls back to the offline
   grounded report card.
5. Click **Deploy**. First build takes a few minutes (installs
   `requirements.txt`). Your app gets a permanent URL like
   `https://agefold.streamlit.app`.

## After deploy

- **Updating the app:** push to `main` and Community Cloud auto-redeploys.
- **First-load latency:** the first time anyone opens a given protein, its
  AlphaFold structure is fetched from EBI and its UniProt/InterPro/ClinVar
  context is pulled live; both are cached afterward.
- **Resource limit:** Community Cloud gives ~1 GB RAM. The parquets load fine;
  if you hit limits later, the graduation path is Docker + a VM or HF Spaces.
- **Rotating the key:** change it in the app's **Settings → Secrets** on
  share.streamlit.io — no redeploy or code change needed.

## Running locally (unchanged)

```bash
cd ~/Documents/lipms-viewer
./run_app.sh          # http://localhost:8501
```
Local runs still use the OneDrive AlphaFold PDBs and your existing
`enrichment_cache.sqlite` when present — the on-demand fetch only kicks in for
proteins whose local file is missing.
