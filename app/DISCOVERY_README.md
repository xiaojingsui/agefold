# Discovery dashboard — ranking worm proteins by aging × disease relevance

The **Discover** tab surfaces which of the 6,823 measured worm proteins are
worth looking at first: proteins whose aging conformational signal overlaps
with pathogenic variants in the human ortholog. That's the paper's central
hypothesis made searchable.

## Scoring formula

```
composite = w_hot × aging_hotspot_max
          + w_load × aging_load_per_aa
          + w_ml × ml_predicted_high_frac
          + alignment_confidence × (
                w_var × variant_load_per_aa
              + w_ovl × aging_variant_overlap
            )
```

| Component | What it measures |
|---|---|
| `aging_hotspot_max` | Max \|log2FC\| among significant aging peptides (day6 + day9, adj-p < 0.05, \|log2FC\| > 1). "At least one strong hotspot." |
| `aging_load_per_aa` | (Number of significant aging peptides × mean \|log2FC\|) / seq_length. "Whole-protein hit, not just one spike." |
| `ml_predicted_high_frac` | Fraction of residues with predicted p(destabilized) > 0.5. Corroborates the observation. |
| `variant_load_per_aa` | Pathogenic ClinVar variants (mapped from top human ortholog) / seq_length. "Human disease relevance." |
| `aging_variant_overlap` | Count of pathogenic variants that fall inside a significant aging residue. **The actually novel signal.** |
| `alignment_confidence` | Worm↔human sequence identity, clipped: <30% → 0, ≥60% → 1, linear in between. Gates the variant component. |

## Default weights

```python
w_hot = 1.0   # peak intensity
w_load = 5.0  # whole-protein load
w_ml = 2.0    # ML corroboration
w_var = 0.5   # variant density
w_ovl = 3.0   # overlap bonus — the actually novel signal
```

Chosen so the top-20 is dominated by known disease genes (SOD1, PTEN, VCP, ACTG1, ribosomal proteins) — a face-validity check. Tune from the app's advanced expander to explore alternative rankings.

## Top-20 hits (default weights)

| Rank | Worm gene | UniProt | Human disease | Composite |
|-----:|-----------|---------|---------------|----------:|
| 1 | *sod-1* | P34697 | Amyotrophic lateral sclerosis type 1 (ALS1) | 149.1 |
| 2 | *fum-1* | O17214 | Hereditary cancer-predisposing syndrome (HLRCC) | 107.6 |
| 3 | *his-1* | P62784 | Tessadori-Van Haaften neurodevelopmental syndrome | 56.1 |
| 4 | *mthf-1* | Q17693 | Homocystinuria | 55.3 |
| 5 | *daf-18* | G5EE01 | PTEN hamartoma tumor syndrome | 55.0 |
| 6 | F02A9.10 | P34385 | 3-methylcrotonyl-CoA carboxylase deficiency | 54.1 |
| 7 | *cdc-48.2* | P54812 | Frontotemporal dementia / ALS (VCP) | 53.3 |
| 8 | *hex-1* | Q22492 | Tay-Sachs disease | 51.9 |
| 9 | *prp-31* | Q9N592 | Retinitis pigmentosa 11 | 45.7 |
| 10 | *rps-19* | O18650 | Diamond-Blackfan anemia | 45.5 |

Every one is a known disease gene with independent OMIM annotation.

## Filter semantics

- **Min alignment identity** — floor for the worm↔human alignment. 0.3 is the recommended default; below that, residue-level correspondence isn't reliable and the variant component is heavily discounted anyway.
- **Min composite score** — cuts the tail. The natural break in the distribution is around 12 (top-50 threshold), where the score histogram cleanly separates.
- **Require ≥1 pathogenic variant mapped** — restricts to proteins whose human ortholog has ClinVar Pathogenic / Likely pathogenic entries. Excludes proteins with only VUS (like HSPA8).
- **Require aging-variant overlap ≥1** — the stricter cut: at least one pathogenic variant falls inside a significant aging residue. Answers "does this protein age in a place that matters for human disease?"
- **Top N** — the table cap. Default 50.

## Known limitations

- **Composite is a ranker, not a p-value.** Nothing here is a statistical test. Use the ranking to prioritize follow-up; use targeted enrichment tests (hypergeometric, permutation) to make quantitative claims about a subset.
- **Top-1 ortholog only.** Where a worm gene maps to multiple human paralogs (actins, tubulins, myosins), the score reflects the top ortholog. All-paralog pooling is a future extension.
- **"Top disease" is per-residue mode, not per-protein.** The disease shown is the most common `top_disease` across variant-carrying residues of that protein. A protein with variants across several disorders (e.g. HSPB1 with two ClinVar-associated CMTs) will show only one.
- **Alignment identity is the weakest link.** At the 30-60% band (majority of the dataset), the variant component is honest but the residue-level correspondence gets qualitative. The pipeline documents this at the residue level via `alignment_identity`; use it.
- **ML component tends to zero for most proteins.** The v2 predictor was trained conservatively; only the strongest aging-signal proteins have residues at p > 0.5. That's why removing the ML component doesn't much change the ranking — the top hits are driven by direct observation × variants, not model belief.
- **Somatic variants are included.** UniProt Variation API returns both germline and somatic ClinVar entries; filter downstream if you need germline-only.

## Reproducing / recomputing

```python
import discovery
df = discovery.compute_scores()                     # defaults
df = discovery.compute_scores({"w_ovl": 10.0})      # emphasize overlap only
```

Output columns: `uniprot_id, gene_symbol, seq_length, composite, aging_hotspot_max, aging_load, aging_load_per_aa, ml_predicted_high_frac, variant_load, variant_load_per_aa, n_res_pathogenic, alignment_identity, aging_variant_overlap, alignment_confidence, top_disease, n_sig_aging_peptides`.

Saved to `app/data/discovery_scores.parquet`; the Streamlit tab reads it via `load_discovery()` and caches it. Rerun the module (or delete the parquet) if the underlying peptide, variant, or ortholog data changes.
