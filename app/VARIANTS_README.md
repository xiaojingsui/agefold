# Disease variant module — methodology and limitations

Maps pathogenic ClinVar variants from human orthologs onto worm residues via pairwise alignment, giving each measured worm protein a residue-level "disease variant density" track alongside its LiP-MS observations and ML predictions.

## Pipeline summary

1. **Ortholog resolution** — For each measured worm protein with a human ortholog (WormBase SimpleMine), take the top ortholog by number of supporting inference methods. Resolve HGNC → human UniProt via UniProt REST search.
2. **Variant fetch** — For each human UniProt, hit the **UniProt Variation API** (`https://www.ebi.ac.uk/proteins/api/variation/{uid}`), which returns ClinVar entries pre-mapped to protein position with clinical significance labels.
3. **Classification** — Variants tagged as "Pathogenic" or "Likely pathogenic" (plus a small number of "Disease"/"Disease causing"/"Association" labels) are marked pathogenic; everything else (Variant of Uncertain Significance, Likely benign, Benign) is not.
4. **Alignment** — Global pairwise alignment of worm sequence (from the local AlphaFold PDB) to the human UniProt sequence using Biopython's `PairwiseAligner` with BLOSUM62 (gap_open=-10, gap_extend=-0.5).
5. **Mapping** — Every human residue in an aligned block gets mapped to its worm counterpart; variants at unmapped positions (in a gap) are dropped.
6. **Per-residue density** — Aggregate to `(worm_uniprot, residue) → (n_variants_all, n_variants_pathogenic, top_disease)`.

## Dataset in numbers

| Metric | Value |
|---|---|
| Measured proteins with any human ortholog | 3,891 |
| Unique human UniProt IDs fetched | 3,420 |
| Human UniProts with any ClinVar variant | 3,400 |
| Human UniProts with ≥1 pathogenic variant | 1,289 |
| Total variant→worm-residue mappings | 2,027,412 |
| Pathogenic mappings | 63,251 |
| Median alignment identity (worm↔human) | 0.41 |
| Alignments passing 30% identity | 3,516 / 3,891 (90%) |
| Worm proteins with ≥1 mapped variant | 3,859 (56.6% of measured) |
| Worm proteins with ≥1 pathogenic mapped | 1,400 (20.5% of measured) |

## Confidence — read this before quoting numbers

Alignment identity is stored per-protein in `variant_alignment_stats.parquet` and per-residue-density row in `variant_density.parquet` (via `alignment_identity`). Treat these tiers as separate:

- **identity ≥ 70 %** — high-confidence 1-to-1 mapping, most residues align. Position of a human variant on the worm sequence is reliable. Actin, tubulin, HSP70/HSP90, ribosomal proteins fall here.
- **30–70 %** — coordinate mapping is qualitative. A pathogenic hotspot in the human ortholog "corresponds to roughly this region" of the worm protein; individual residue-level claims are speculative. The majority of the dataset (~median 41 %) sits here.
- **< 30 %** — treat as "there is a related human protein with disease variants, but you cannot claim residue-level correspondence." The pipeline still stores the mapping but you should not use it for causal claims.

## Known limitations

- **VUS excluded by default.** Most C. elegans orthologs are essential and highly conserved, so ClinVar entries on their human counterparts tend to be dominated by VUS (variants of uncertain significance). HSP-1 → HSPA8 is the canonical example: 588 variants in ClinVar, 0 pathogenic. That's a real biological fact — HSPA8 has no OMIM disease associations — not a mapping failure.
- **Top-1 ortholog only.** For worm genes with multiple human paralogs (e.g. worm ACT-1/2/3/4/5 all map to actin paralogs), we take the top one by SimpleMine method count. All-paralog pooling is a natural extension.
- **Alignment gaps drop residues.** Variants at insertion/deletion positions in the human sequence relative to the worm cannot be mapped and are silently dropped from the density; this affects a small fraction of variants (median alignment fills ~92% of residues, based on the identity-vs-mapped-count scatter).
- **No structural rescue.** Where the sequence-level alignment fails, we don't try to rescue with structural superposition. That would help for the low-identity end of the distribution.
- **API cache is authoritative.** All fetched variant and HGNC payloads sit in `app/data/enrichment_cache.sqlite` (source keys `uniprot_variation`, `hgnc`). Rebuilding the density requires re-running the alignment pass, not re-fetching from the API.
- **Somatic variants included by default.** UniProt Variation API returns both germline and somatic entries; the pipeline does not filter by `somaticStatus`. Filter downstream if you need germline-only.
- **TF ortholog caveat carries over.** SimpleMine can misassign orthologs for fast-evolving TFs (SKN-1 → TCP1 example in MODEL_CARD §7). The variant module inherits any ortholog error the upstream table makes.

## Files

```
app/
├── variants.py                       ← fetch + align + map helpers
├── data/
│   ├── variant_density.parquet       ← per-(worm_uniprot, residue) counts + top disease
│   └── enrichment_cache.sqlite       ← backing cache for UniProt Variation
data/                                 ← build-only, not shipped in app
├── human_variants_by_worm.parquet    ← raw 2M-row mapping (every variant × every mapped residue)
├── variant_alignment_stats.parquet   ← per-protein alignment identity + n_mapped
├── worm_to_human_resolution.parquet  ← worm_uniprot → top HGNC → human UniProt
└── variant_probe.json                ← small probe on HSPA8 + ACTB
```

## Where this feeds the app

- **Streamlit tracks** — `streamlit_app.py` renders a fourth track under the ML prediction track, colored by pathogenic variant density (Reds cmap, 95th-percentile clip).
- **Chat retrieval** — `retrieval.py::get_protein_context()` returns a `variant_overlap` section with per-aging-region variant counts and top diseases. The system prompt (rule 8) requires the LLM to cite these as `[S7]` and distinguish them from LiP-MS observations.
- **Right-column card** — "Disease variants" panel shows top diseases and the 5 most variant-dense residues on the current protein.
