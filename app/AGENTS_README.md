# Multi-agent analysis framework

The **Analysis workspace** (the app's default tab) turns the viewer into a
conversational, multi-agent structural-proteomics assistant. You ask a question
in plain language; a coordinator routes it to the relevant specialists; each
reasons over the grounded data; and their answers are synthesized into one cited
response — with a visible trace of which agents contributed.

This layer sits **upstream** of everything else. The manual viewer (Overview),
single-protein chat, Discover dashboard, ML track, variant track, structural
context, and the public/researcher dual-track are all preserved and reachable.

---

## Architecture

```
                         ┌──────────────────────────┐
   your question  ─────▶ │   Coordinator / Router   │   agents.route()
                         └──────────────────────────┘
                                     │ selects 1–3 specialists
             ┌───────────────┬───────┴───────┬────────────────┐
             ▼               ▼               ▼                ▼
     🧬 Structural   📈 Aging & Stress  🩺 Disease &     📊 Data
        Biologist        Analyst         Genetics         Analyst
      [S8·S2·S3]          [S5]           [S4·S7]        [S6 + dataset]
             └───────────────┴───────┬───────┴────────────────┘
                                     ▼
                         ┌──────────────────────────┐
                         │        Synthesis         │  one cited answer
                         └──────────────────────────┘   + "agents consulted" trace
```

Three files:

| File | Role |
|------|------|
| `agents.py` | Specialist **registry** (roles, prompts, focus, sources) + keyword **router** (`route`, `route_explain`, `is_dataset_level`). Pure data/logic — no network, no Streamlit. |
| `analysis_tools.py` | **Dataset-level** functions over the precomputed parquets (cross-protein / cross-condition). |
| `agent_orchestrator.py` | `ask_agents()` — route → run specialists in parallel → synthesize → stream. Reuses `chat_backend` key/client and `prompts.format_context_prompt`. |

---

## The specialists

Each specialist sees the **same grounded context** (`retrieval.get_protein_context`,
sources `[S1]`–`[S8]`) but is told to stay in its lane and centers on its own
sections.

| Agent | Scope | Focus sections | Sources |
|-------|-------|----------------|---------|
| 🧬 **Structural Biologist** | Fold stability, buried/exposed, secondary structure, pLDDT, distance to functional sites | `structural_context`, `domains`, `aging_signal`, `function` | `[S8]` `[S2]` `[S3]` `[S1]` |
| 📈 **Aging & Stress Analyst** | The 9-condition comparison — aging-specific vs shared with heat-shock / polyQ / ts-mutant | `aging_signal`, `lipms`, `protein_annotations` | `[S5]` |
| 🩺 **Disease & Genetics** | Human orthologs + ClinVar variant mapping, with homology caveats | `orthologs`, `variant_overlap`, `identity` | `[S4]` `[S7]` `[S1]` |
| 📊 **Data Analyst** | ML destabilization scores + **dataset-level** queries (cross-protein/condition) | `ml`, `aging_signal`, + dataset tools | `[S6]` `[S5]` |

**Shared grounding rules** (all specialists): answer only from the provided
context; cite every claim with a bracketed source; keep the evidence-type
discipline — measured LiP-MS `[S5]` vs ML prediction `[S6]` vs homology inference
`[S7]` vs AlphaFold-predicted structure `[S8]`; prefer specific residues/regions;
be concise (one voice on a panel).

---

## The router

`agents.route(question)` scores the question against per-specialist keyword sets
and returns the **1–3** highest-scoring specialists (default cap 3). Rules:

- **Dataset-level** questions ("across all proteins", "which proteins", "top…",
  "how many", "enrich") add the **Data Analyst** with a boost.
- **No keyword hit** → a sensible default panel (Structural Biologist + Aging &
  Stress Analyst) for an open-ended protein question.
- `route_explain()` returns the per-specialist scores and the dataset-level flag
  — this is what the UI shows as the routing trace.

Routing is deterministic and instant (no LLM call), so it's cheap and testable.

---

## Dataset-level analysis (`analysis_tools.py`)

The Data Analyst can answer **cross-cutting** questions the single-protein view
can't, all from the precomputed parquets:

| Function | Answers |
|----------|---------|
| `compare_conditions(uid)` | Peak \|log2FC\| + significant-peptide counts per condition, grouped aging / stress / polyQ / ts-mutant, with an aging-specific-vs-shared verdict. |
| `top_movers(condition, n)` | Proteins with the largest peak \|log2FC\| in a condition. |
| `aging_specific_hits(n, min_ratio)` | Proteins whose aging peak is high **and** clearly exceeds their stress/polyQ/ts peak. |
| `structural_enrichment_summary()` | Whole-dataset: are aging hotspots enriched in buried / near-site regions? (the 1.31× near-site finding). |
| `discovery_top(n)` | Top composite-scored disease-linked candidates. |

Each returns `{available, table, summary}`; the orchestrator injects the
`summary` lines into the Data Analyst's context as a `DATASET ANALYSIS` block.

---

## Cost & latency

Per question (online mode, Sonnet):

- **Router:** free (no LLM).
- **Specialists:** 1 call each, in **parallel** — so wall-time ≈ one specialist
  call (~4–8 s), not the sum. Routing to 1–3 keeps this bounded.
- **Synthesis:** 1 streamed call (~4–8 s).
- **Total:** typically **2–4 LLM calls**, ~10–20 s, roughly **$0.03–0.12** per
  turn depending on how many specialists fire and answer length. Specialists are
  capped at 900 output tokens, synthesis at 1400.

Compared to the single-protein chat (1 call), the multi-agent path costs more but
gives specialist depth, cross-condition/dataset reasoning, and a transparent
trace.

---

## Offline mode

With no Anthropic API key, `ask_agents` still routes (showing which specialists
*would* answer) and renders the grounded **report card** from the context — so
the workspace degrades gracefully instead of erroring. Set `ANTHROPIC_API_KEY`
or `.streamlit/secrets.toml` to enable the live panel.

---

## Adding a specialist

1. Add an entry to `AGENTS` in `agents.py`: `name`, `icon`, `role`, `focus`
   (context sections), `sources` (citation codes), and a `system_prompt` that
   starts from `_SHARED_RULES` and adds the new lane.
2. Add the id to `SPECIALIST_ORDER`.
3. Add keyword signals to `_ROUTER_KEYWORDS[<id>]`.
4. (Optional) If it needs dataset-level data, add a function to
   `analysis_tools.py` and wire it into `_gather_dataset_results`.

No orchestrator or UI change is needed — the workspace reads the registry.

---

## Limitations

- **Grounding is only as good as the context.** Specialists answer from
  `get_protein_context`; gaps in UniProt/InterPro/variant data propagate.
- **Router is keyword-based.** It's fast and predictable but can miss an
  unusually phrased question; it errs toward including an extra specialist rather
  than too few.
- **Homology caveats stand.** Disease-variant claims are ClinVar entries on the
  human ortholog mapped by alignment — inferences, not worm observations.
- **Structure is AlphaFold-predicted**, monomer geometry — a residue buried at a
  real oligomeric interface can read as exposed.
- **Dataset tools are precomputed lookups**, not live re-analysis; they reflect
  the parquets shipped in `app/data/`.
