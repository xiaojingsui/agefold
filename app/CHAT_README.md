# Chat interface — setup and design notes

The **Ask about this protein** tab in the viewer lets users ask free-form questions about the currently loaded protein. Every claim in the answer is cited to one of six structured sources.

## Quick start — enabling chat

The chat feature is optional and disabled by default. Two ways to enable it:

### Option 1 — environment variable (recommended for local dev)

```bash
export ANTHROPIC_API_KEY=sk-ant-…
streamlit run app/streamlit_app.py
```

### Option 2 — Streamlit secrets (recommended for deployment)

Create `.streamlit/secrets.toml` in the app root:

```toml
anthropic_api_key = "sk-ant-…"
```

Get a key at [console.anthropic.com](https://console.anthropic.com/) → API Keys. **Never commit the key** — `.streamlit/secrets.toml` is in `.gitignore` by convention.

## Cost per turn

The chat calls **claude-sonnet-4-5** by default. Per turn:
- **Input:** ~3,000–5,000 tokens (structured context is ~4 kB for a typical protein)
- **Output:** ~800–1,200 tokens for a substantive question
- **Approximate cost:** $0.02–$0.04 per question at current Sonnet pricing (Sep 2024)

Long DAF-16-style questions on proteins with a lot of context can push input tokens closer to 8k. See [Anthropic pricing](https://www.anthropic.com/pricing).

## What happens without a key

The tab shows an **offline report card** — a rendered markdown summary of the structured context: identity, function, domains, aging hotspots, ML top-predicted residues, and orthologs/diseases. It's not conversational but it's still informative. Users can read it to get 80% of what the chat would answer for the "what is this protein and what changes with age" archetype question.

## How retrieval works

`app/retrieval.py::get_protein_context(uid)` assembles a nested dict from these sources:

| Section | Source | File |
|---------|--------|------|
| identity | UniProt REST + LiP-MS Table S3 | `protein_features.parquet` + `enrichment_cache.sqlite` |
| function | UniProt REST (cached in SQLite) | `enrichment_cache.sqlite` |
| domains | InterPro (bulk C. elegans dump) | `interpro_domains.parquet` |
| lipms | Table S1 peptide-level | `lipms_peptide_unified.parquet` |
| aging_signal | Computed from lipms — significant peptides in day6/day9 grouped into ≤20-aa contiguous regions | (derived) |
| ml | LightGBM predictor v2 (baseline + ESM-2 PCs) | `residue_predictions_v2.parquet` |
| orthologs | WormBase SimpleMine (parsed at build time) | `orthologs_and_disease.parquet` |
| sources | Hard-coded URL keys S1-S6 | (in code) |

Every prompt sent to the LLM includes the full sources block and the system prompt enforces that no claim can be made without a bracketed citation.

## Prompt design

`app/prompts.py` implements two things:

1. **`SYSTEM_PROMPT`** — 7 rules covering citation, grounding, evidence-type distinction (observed LiP-MS vs ML predictions), specificity, structural interpretation, brevity, and fallback for general questions.
2. **`format_context_prompt(ctx, question)`** — Jinja-style template that flattens the retrieval dict into a ~4 kB structured message.

## Known limitations

- **UniProt function is missing for some proteins** (e.g. HSP-1/P09446). The chat handles this gracefully by falling back to GO annotations and InterPro domain reasoning, but the answers are less complete when this happens.
- **DAF-16 and other transcription factors** have very sparse LiP-MS coverage (2 peptides for DAF-16, 0 for HSF-1, SKN-1, PHA-4). The dataset targets abundant proteins. The chat explicitly acknowledges when it can't distinguish "no aging change" from "insufficient coverage."
- **SimpleMine ortholog mapping can be wrong for TFs** (see MODEL_CARD.md §7 — SKN-1 maps to TCP1 in SimpleMine, which is a known WormBase artifact for fast-evolving TFs). Users should sanity-check TF ortholog claims against primary literature.
- **Disease associations are inferred "by orthology"** — these come from human data mapped through the worm→human ortholog table. They are hypotheses, not evidence that the C. elegans gene is disease-associated.
- **ML predictions are prioritization scores, not measurements.** The prompt makes the model distinguish these explicitly, but users should read `MODEL_CARD.md` for the failure modes (especially on low-abundance proteins like HSP-1 and ACT-2).
- **No live literature search.** The chat does not retrieve PubMed abstracts, so novel or recent findings will be missing.
- **No structural analysis at atomic resolution.** The chat can say "residues 454-500 are in domain X" but it cannot say "these residues are on the domain interface" — that would need direct AlphaFold parsing.

## Example transcripts

See `chat_examples.md` in this folder — 5 proteins × 1-2 questions each, ~800-1,700 output tokens per answer. In the 8-turn evaluation set, every response ended cleanly (`stop_reason=end_turn`) and cited at least one source from the [S1]-[S6] block (average 3.4 distinct sources per answer). Grounding is not perfect — the model can still make interval-assignment errors when a region straddles two adjacent domains — but the retrieval layer pre-computes domain overlaps for each aging region, so region→domain mapping is done deterministically upstream rather than left to the LLM.

## Files this feature depends on

```
app/
├── retrieval.py           ← builds the structured context
├── prompts.py             ← system prompt + template
├── chat_backend.py        ← Anthropic SDK client with offline fallback
├── streamlit_app.py       ← the chat tab UI
├── chat_examples.md       ← 5-protein QA transcripts
└── data/
    ├── protein_features.parquet
    ├── lipms_peptide_unified.parquet
    ├── interpro_domains.parquet
    ├── residue_predictions_v2.parquet
    ├── orthologs_and_disease.parquet
    ├── enrichment_cache.sqlite
    └── alphafold_index.parquet
```
