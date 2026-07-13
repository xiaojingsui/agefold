# AgeFold — a conversational AI to decode the biology of aging

An interactive resource for the *C. elegans* TMT-LiP-MS aging dataset. This
repository builds the **first component** of the larger project: a per-protein
viewer that paints per-residue LiP-MS conformational signal onto AlphaFold
structures, alongside domain, feature, and homology context.

The ML feature-prediction module and the conversational chat module both
plug into the same backend — see [`app/DATA_SCHEMA.md`](app/DATA_SCHEMA.md)
for the contract they read against.

---

## Analysis workspace — multi-agent framework (default surface)

The app opens on a conversational **🧪 Analysis** workspace. Ask a
structural-proteomics question in plain language; a **coordinator** routes it to
the relevant specialists — **🧬 Structural Biologist**, **📈 Aging & Stress
Analyst**, **🩺 Disease & Genetics**, **📊 Data Analyst** — which each reason over
the grounded data in their lane, and their answers are **synthesized into one
cited response** with a visible "agents consulted" trace. The Data Analyst can
also answer **dataset-level** questions (aging-specific hits, top movers,
structural enrichment across all proteins), not just single-protein ones.

This framework is the new **upstream entry point**; every existing surface — the
manual viewer (Overview/Discover), ML track, variant track, structural-context
card, and the researcher/public dual-track — is preserved and reachable. See
[`app/AGENTS_README.md`](app/AGENTS_README.md) for the architecture, each agent's
scope and sources, cost/latency, and how to add a specialist.

---

## Two doors: Researcher and Explore-aging

A **View** toggle at the top of the sidebar switches the whole app between two
audiences reading the *same* data:

- **🔬 Researcher** — the full technical viewer: 9-condition per-residue heatmap,
  ML destabilization track, mapped disease variants, InterPro domains, and chat
  with `[S1]`–`[S7]` citations.
- **🌍 Explore aging** — a plain-language tour for the curious public: a landing
  page with a how-to-read guide and a curated gallery of 7 featured proteins,
  each re-narrated as an accessible aging story (what it does, where it changes
  with age, human-health connection, why it matters), plus a jargon-free chat.

The public track is **re-narration, not dumbing-down**: identical measurements
and retrieval context, plain language on top, every claim verified against the
full dataset. See [`app/PUBLIC_README.md`](app/PUBLIC_README.md) for the design
philosophy, curation criteria, and accuracy safeguards.

Implementation: `public_view.py` (renderer) + `public_content.py` (curated,
fact-checked stories) + `SYSTEM_PROMPT_PUBLIC` in `prompts.py`. The public branch
runs then `st.stop()`s, so all researcher code is untouched.

---

## What's in the box

- `app/streamlit_app.py` — the app (protein page)
- `app/enrichment.py` — cached UniProt / InterPro / HGNC lookups
- `app/data/*.parquet` — the six unified data tables
- `app/DATA_SCHEMA.md` — full column dictionary and recipes
- `app/README.md` — how to run the app

## Run it

```bash
conda activate lipms-viewer
cd app
streamlit run streamlit_app.py
```

## What the protein page shows

For any UniProt ID or gene symbol you type:

1. **Interactive 3D AlphaFold cartoon**, painted per-residue by log₂FC on a
   zero-centered diverging colormap (blue = destabilized, red = stabilized;
   unmeasured residues grey).
2. **Nine-condition residue heatmap** (day 6 / day 9 aging, heat shock,
   Q35 / Q40 polyQ, myosin-ts 15/25°C, paramyosin-ts 15/25°C), current
   condition boxed.
3. **InterPro domain track** aligned to the heatmap x-axis.
4. **Protein card**: UniProt name and function (live), subcellular location,
   sequence length, MW, log₂(emPAI), network degree centrality, longevity
   flag, ribosome pausing, turnover, ubiquitination, HSP90 / TRiC client
   status, tissue specificity, GO BP/MF/CC, keywords.
5. **Per-peptide table** for the current condition, sorted by |log₂FC|.
6. **Worm → human ortholog panel**, with UniProt/HGNC links and inferred
   disease associations from WormBase SimpleMine.

## Dataset in numbers

- 6,823 proteins with LiP-MS measurements across 9 conditions
- 530,947 peptide × condition measurements
- 94.7% of measured proteins have a local AlphaFold structure
- 79.7% have InterPro domain annotations
- 58% have a curated worm → human ortholog
- 23% have a disease association via that ortholog

## Roadmap

| Module              | Status | Depends on |
|---------------------|--------|-----------|
| Protein page (v0)   | ✅ shipped | — |
| ML residue predictor v2 (LightGBM + ESM-2 PCs) | ✅ shipped — AUC 0.699, AP 0.331 | `lipms_per_residue_long.parquet` + `protein_features.parquet` |
| Conversational chat with cited answers | ✅ shipped — see `app/CHAT_README.md` | Anthropic API key (env or `.streamlit/secrets.toml`) |
| ClinVar variant mapping onto worm structure | ✅ shipped — see `app/VARIANTS_README.md` | UniProt Variation API + pairwise alignment |

## Chat interface

The "Ask about this protein" tab lets users ask free-form questions grounded in
UniProt, InterPro, AlphaFold, WormBase, the LiP-MS dataset, the ML predictor,
and (as of the disease variant module) ClinVar via UniProt Variation. Every
claim in the answer is cited to one of seven structured sources [S1]-[S7],
with the system prompt explicitly requiring the model to distinguish observed
LiP-MS changes from ortholog-mapped disease variants.

Without an API key, the tab shows an **offline report card** — a rendered
markdown summary of the protein's structured context (still useful; no
conversation, but no hallucination either).

See [`app/CHAT_README.md`](app/CHAT_README.md) for key setup, prompt design,
retrieval sources, and known limitations. See
[`app/chat_examples.md`](app/chat_examples.md) for example transcripts on
VIT-6, HSP-1, ACT-2, UNC-54, and DAF-16.

## Discovery dashboard

`app/discovery.py` computes a per-protein composite score from four signals:
peak aging \|log2FC\|, aging load per residue, ML confidence, and variant overlap
with the human ortholog. The **🔍 Discover** tab exposes this as a sortable
top-hits table with filters (alignment identity, min score, require variants,
require aging-variant overlap).

**Default top-20 is dominated by known disease genes** — SOD-1 → ALS1,
daf-18 → PTEN, cdc-48.2 → VCP-FTD-ALS, prp-31 → retinitis pigmentosa,
rps-19 → Diamond-Blackfan, ACT-2 → hearing loss. That's the paper's central
hypothesis made face-valid in one screen. See
[`app/DISCOVERY_README.md`](app/DISCOVERY_README.md) for the scoring formula,
weight defaults, and known limitations.

## Structural context module

`app/structure_features.py` reads the local AlphaFold model for every measured
protein and computes, per residue: solvent accessibility (buried core vs exposed
surface), three-state secondary structure (helix/sheet/loop), pLDDT confidence,
and — where UniProt annotates a functional site — the CA–CA distance to the
nearest active/binding/metal/disulfide site. This turns *"the structure changes
here"* into *"this buried, functionally coupled region changes — the kind of
change most likely to affect fold stability."*

Shown as a **Structural context** card in the researcher panel (per aging region:
buried/surface, secondary structure, confidence dot, distance to site, plain
interpretation) and threaded into the chat as citation source **`[S8]`**.

**Coverage:** 6,462 proteins / 3.26M residues annotated (0 errors); 1,409
proteins carry UniProt functional sites (distance-to-site computed for 1,341 with
a local structure). Whole-dataset finding: aging hotspots
are structurally similar to background overall (~38% buried both) but **1.31×
enriched near functional sites**. Click-to-couple lets you select a residue range
and see it highlighted on both the 3D structure and the heatmap at once. See
[`app/STRUCTURE_README.md`](app/STRUCTURE_README.md) for methods, thresholds,
and limitations.

## Disease variant module

`app/variants.py` fetches ClinVar variants for each protein's top human ortholog
via the UniProt Variation API, pairwise-aligns worm↔human sequences with BLOSUM62,
and maps every variant onto its worm-residue counterpart. The result is a per-residue
"pathogenic variant density" track shown under the ML prediction track and a
"Disease variants" panel in the right column, plus a new `[S7]` citation source for
the chat.

**Coverage:** 3,859 / 6,823 measured proteins have at least one mapped variant
(56.6%), 1,400 have ≥1 pathogenic variant (20.5%). Median alignment identity is
0.41. See [`app/VARIANTS_README.md`](app/VARIANTS_README.md) for the full
methodology and confidence tiers.

## ML predictor

`app/predictor.py` scores every observed residue with a probability of being
"significantly destabilized in aging" (|log2FC|>1 AND adj-p<0.05 in day6 or
day9). Best model uses baseline structural/positional features plus ESM-2 35M
per-residue embeddings (compressed to 32 PCs). See
[`app/MODEL_CARD.md`](app/MODEL_CARD.md) for the ablation table, calibration
notes, and seven documented failure modes.

## Build summary

Four figures produced during the build (also saved as artifacts):

- `data/qc_coverage.png` — dataset coverage across 9 conditions
- `data/vit6_residue_track_qc.png` — residue-track QC on VIT-6
- `data/app_preview_vit6_full.png` — full app protein page preview (VIT-6, day 6)
- `data/app_preview_orthologs.png` — ortholog panel preview (HSP-1, SKN-1)
- `data/build_summary_gallery.png` — the four combined
