# Public "Explore aging" track — design & maintenance

The app has two doors into the **same data**, chosen with the **View** toggle at
the top of the sidebar:

- **🔬 Researcher** — the full technical viewer (9-condition heatmap, ML
  destabilization track, mapped disease variants, InterPro domains, cited chat).
- **🌍 Explore aging** — a plain-language tour for the curious public.

This document explains the public track: its philosophy, how the featured set is
curated, the accuracy safeguards, and how to extend it.

---

## Design philosophy: re-narration, not dumbing-down

The public track is **not** a stripped-down researcher view. It is the same
measurements, the same 3D structures, and the same underlying retrieval context —
*re-narrated* for someone with no biology background. Nothing is invented and
nothing is exaggerated to make a better story. Where the data is uncertain, the
public copy says so in plain words ("this is based on the human cousin of the
gene, not measured in people directly").

The guiding test for every sentence of public copy: **would a domain expert
reading it agree it is accurate, even though it uses no jargon?**

---

## What the public track contains

| Piece | File | Notes |
|-------|------|-------|
| Mode toggle | `streamlit_app.py` | `st.session_state["view_mode"]`; public branch runs `public_view.render_public()` then `st.stop()`, leaving researcher code untouched |
| Landing page | `public_view.py` `_render_landing` | Headline, intro, how-to-read diagram, why-*C. elegans*, featured gallery, glossary |
| How-to-read diagram | `data/public_howto_diagram.png` | Colour legend + "spin the structure" schematic |
| Featured stories | `public_content.py` `FEATURED` | 7 curated proteins (pure data + text) |
| Glossary | `public_content.py` `GLOSSARY` | 10 plain-language definitions |
| Protein view | `public_view.py` `_render_protein` | 3D structure, "what it does", "where it changes", human-health, "why it matters", researcher escape-hatch |
| Public chat | `public_view.py` `_render_public_chat` + `prompts.py` `SYSTEM_PROMPT_PUBLIC` | Same retrieval context; lay-audience system prompt; no `[S#]` codes shown |

---

## Curation criteria for featured proteins

A protein earns a spot in `FEATURED` only if **all** of these hold:

1. **Real aging signal in this dataset.** It must have significant conformational
   change in day6/day9 (`|log2FC| > 1`, `adj-p < 0.05`). We dropped `daf-16`
   (the famous longevity gene) precisely because it shows *no* aging signal here —
   there is no honest "changes with age" story to tell about it.
2. **A compelling, true hook.** Either a clear human-disease link (sod-1→ALS,
   daf-18→PTEN cancer, cdc-48.2→VCP dementia, act-2→hearing loss) or an
   intrinsically vivid biology story (vit-6 yolk build-up, unc-54 muscle motor,
   hsp-1 the repair crew).
3. **Coverage of different aging themes.** The set spans oxidative damage,
   protein quality control, muscle decline, the cytoskeleton, growth/longevity
   signalling, and reproduction — so the gallery teaches the breadth of aging
   biology, not one narrow slice.

Human-health sections appear **only** when there is a real worm→human ortholog
with pathogenic ClinVar variants. Pure-biology stories (vit-6, unc-54, hsp-1)
omit that section rather than stretch for a weak link.

---

## Accuracy safeguards

- **Every factual claim in `public_content.py` was cross-checked** against
  `retrieval.get_protein_context()` and the UniProt cache before shipping.
- **Superlatives are checked against the full dataset**, not a hand-picked
  sample. Example: vit-6's copy says it is "among the very largest" age-related
  changes — verified as rank 30 of 5,252 proteins by peak aging `|log2FC|`, so we
  do **not** claim it is the single largest (perm-4 leads at `|log2FC|`≈9.6).
- **Spatial claims are computed.** The sod-1 story says human ALS mutations
  "cluster on the same aging region (82–137)" — verified: 29 of 74 pathogenic
  worm residues fall inside 82–137 (`variant_density.parquet`, alignment identity
  0.59).
- **The public chat uses the identical retrieval context** as the researcher
  chat. `SYSTEM_PROMPT_PUBLIC` changes only the *language and tone* — the grounding
  rules ("never invent protein-specific facts", "say when the data doesn't cover
  something") are preserved. Live-tested: zero jargon leaks across function,
  disease, and where-it-changes questions on sod-1 and vit-6.
- **Inference is always flagged.** Disease connections are described as coming
  from "the human cousin of this gene" and noted as similarity-based, never as a
  measurement in humans.

---

## How to add a featured protein

1. Confirm it has real aging signal:
   ```python
   import retrieval
   ctx = retrieval.get_protein_context("<UNIPROT_ID>")
   print(ctx["aging_signal"]["top_regions"])  # must be non-empty
   ```
2. Pull its function and any disease link from the same context (and
   `enrichment.fetch_uniprot`). Do **not** write from memory.
3. Add an entry to `FEATURED` in `public_content.py` with: `uid`, `gene`, `icon`,
   `title`, `hook`, `story`, `why_matters`, and `human_health` **only if** there
   is a real ortholog+pathogenic-variant link.
4. If you make a superlative or spatial claim, verify it against the full
   dataset / `variant_density.parquet` first (see Accuracy safeguards).
5. Headless-check both modes still import:
   ```bash
   python -c "import importlib.util,sys,os; sys.path.insert(0,'app'); \
     importlib.util.spec_from_file_location('m','app/streamlit_app.py').loader.exec_module(\
     importlib.util.module_from_spec(importlib.util.spec_from_file_location('m','app/streamlit_app.py')))"
   ```

---

## Known limitations

- Featured stories are hand-written; adding one is a manual, verified step (by
  design — accuracy over volume).
- Non-featured proteins opened in public mode fall back to a cleaned UniProt
  function string, which can be terse or missing for poorly-annotated proteins.
- The public chat needs an Anthropic API key (same as the researcher chat).
  Without one, the protein stories and 3D view still work; only the free-form
  Q&A is disabled.
- Disease links rely on worm→human ortholog alignment; for low-identity
  orthologs the connection is qualitative. The featured set deliberately favours
  high-identity, well-established orthologs.
