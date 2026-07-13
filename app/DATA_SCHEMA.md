# DATA_SCHEMA — LiP-MS Aging Structure Viewer

Schema and file dictionary for the app backend. This is the contract the ML
module and the conversational chat module read against; extend but don't
break existing columns.

---

## Overview

Source dataset: **TMT-LiP-MS in *C. elegans*** across nine conditions,
comparing structural stability at the peptide level. Nine measured
conditions collapse into four biological groups:

| Group      | Conditions                                             | Comparison base       |
|------------|--------------------------------------------------------|-----------------------|
| aging      | `day6`, `day9`                                         | vs WT day 1           |
| stress     | `hs` (heat shock 35 °C)                                | vs WT 20 °C           |
| polyQ      | `q35`, `q40`                                           | vs Q24                |
| ts-mutant  | `myosin_ts_15/25`, `paramyosin_ts_15/25`               | vs WT at same temp    |

**Dataset shape**: 108,070 peptides · 6,823 proteins · ~530,947 measured
(peptide × condition) rows. 94.7 % of measured proteins have a local
AlphaFold structure (v2, UP000001940).

**Coordinate convention**: 1-based UniProt residue numbering throughout.
Every AlphaFold PDB in the bundle numbers CA atoms 1…seq_length with no
offset, verified end-to-end on VIT-6 (1,651 residues).

---

## Parquet files (in `data/`)

### `lipms_peptide_unified.parquet` — 530,947 rows

The **primary long-form measurement table** — one row per (peptide, condition).

| Column                  | Type    | Description |
|-------------------------|---------|-------------|
| `peptide_sequence`      | str     | Peptide, `X.SEQVENCE.Y` (flanking residues around cut site). |
| `protein_id`            | str     | Full "sp\|ACCESSION\|NAME_CAEEL" or "tr\|…" tag. |
| `uniprot_id`            | str     | Bare UniProt accession (`P18948`). **Primary join key.** |
| `gene_symbol`           | str     | WormBase public name (`vit-6`, `hsp-1`). |
| `start_pos`, `end_pos`  | Int64   | Residue span, 1-based inclusive. |
| `peptide_type`          | str     | `full` / `half` / `semi` — how many PK-cut termini. |
| `condition`             | str     | See condition vocabulary below. |
| `log2fc`                | float64 | Average log₂ fold-change (LiP-MS conformational signal). |
| `pval`                  | float64 | Limma p-value. |
| `adj_pval`              | float64 | BH-adjusted p-value. |
| `pk_cut_residue`        | str     | The amino acid at the PK cleavage site. |
| `avg_residue_depth`     | float64 | Mean residue burial from AlphaFold structure. |
| `secondary_structure`   | str     | DSSP-derived (`H`, `E`, `T`, `-`, …). |
| `idr_annotation`        | str     | `idr` if peptide lies in an intrinsically-disordered region. |
| `interpro_id`           | str     | InterPro accession the peptide overlaps (may be null). |
| `interpro_name`         | str     | Human-readable domain name. |

### `lipms_per_residue_agg.parquet` — 3.79 M rows

**Residue-level aggregation for structure painting.** One row per
(protein, condition, residue). Where multiple peptides cover a residue,
`log2fc_max_abs` picks the peptide with the largest |log2FC| (the most
striking conformational signal at that position).

| Column                     | Type   | Description |
|----------------------------|--------|-------------|
| `uniprot_id`, `gene_symbol`, `condition` | str | — |
| `residue`                  | int32  | 1-based position. |
| `log2fc_max_abs`           | float  | log₂FC of the peptide with max \|log2FC\| covering this residue. |
| `log2fc_mean`              | float  | Mean log₂FC across covering peptides. |
| `adj_pval_at_max_abs`      | float  | Adj-p of the picked peptide. |
| `peptide_type_at_max_abs`  | str    | `full`/`half`/`semi` at the picked peptide. |
| `n_peptides`               | int32  | Number of covering peptides for this (protein, condition, residue). |

### `lipms_per_residue_long.parquet` — 5.74 M rows

**Peptide-per-residue explosion for ML.** Same shape as the agg table, but
one row per (peptide, residue) so the peptide identity is preserved. Use
when you need to trace which peptide contributes each residue's signal.

### `protein_features.parquet` — 9,509 rows

Per-protein annotations from manuscript Table S3.

| Column                    | Type | Description |
|---------------------------|------|-------------|
| `uniprot_id`, `protein_id`, `gene_symbol` | str | — |
| `seq_length`              | int   | Amino-acid length. |
| `mol_weight`              | float | Da. |
| `charge`                  | float | Net charge at pH 7. |
| `log2_empai`              | float | log₂(exponentially modified protein abundance index) — expression proxy. |
| `degree_centrality`       | float | Network degree centrality in the C. elegans PPI graph. |
| `hsp90_client`            | str   | If HSP90 client, source (Taipale/Karagöz/…). |
| `tric_client`             | str   | If TRiC client. |
| `tissue_specificity`      | str   | e.g., `hypodermis`, `neuron`, `not tissue-specific`. |
| `longevity`               | str   | `Pro-Longevity` / `Anti-Longevity` / null (from GenAge or manuscript curation). |
| `age_ribosome_pausing`    | str   | Age-dependent ribosome-pausing category. |
| `ubiquitination`          | str   | Ubiquitination status. |
| `age_protein_turnover`    | str   | Age-dependent turnover category. |
| `measured_in_lipms`       | bool  | True if present in the peptide-level table. |

### `interpro_domains.parquet` — 57,586 rows

InterPro domain intervals for 15,794 *C. elegans* proteins.

| Column          | Type   | Description |
|-----------------|--------|-------------|
| `uniprot_id`    | str    | Uppercase. |
| `interpro_id`   | str    | IPR accession. |
| `interpro_name` | str    | Domain name. |
| `domain_type`   | str    | `domain` / `family` / `homologous_superfamily` / `repeat` / `conserved_site` / `active_site` / `binding_site` / `ptm`. |
| `domain_start`, `domain_end` | int | 1-based inclusive. |
| `protein_length`| int    | For sanity-checking. |
| `dc_status`     | str    | InterPro's domain-coverage status. |

### `alphafold_index.parquet` — 19,694 rows

Local AlphaFold-v2 PDB path lookup.

| Column       | Type | Description |
|--------------|------|-------------|
| `uniprot_id` | str  | — |
| `pdb_path`   | str  | Absolute path to a local `AF-<uid>-F1-model_v2.pdb`. |

### `orthologs_and_disease.parquet` — 5,911 rows

WormBase SimpleMine worm→human ortholog mapping with disease associations.

| Column                    | Type | Description |
|---------------------------|------|-------------|
| `worm_uniprot`            | str  | Primary join key. |
| `worm_gene`, `wormbase_id`| str  | — |
| `orthologs`               | str (JSON list) | `[{"hgnc_id": "HGNC:5241", "n_methods": 7, "methods": "OMA;Compara;Panther;..."}]` |
| `n_orthologs`             | int  | Length of the list. |
| `diseases`                | str (JSON list) | `[{"term": "...", "source": "By Orthology"}]` |
| `n_diseases`              | int  | — |
| `concise_description`, `automated_description` | str | WormBase gene descriptions. |

Human UniProt for each HGNC is resolved on demand via the REST API and
cached in the SQLite database below.

### `enrichment_cache.sqlite`

Runtime cache written by `enrichment.py`. Persists across app reboots.

```
enrichment (uniprot_id TEXT, source TEXT, payload TEXT JSON, fetched_at REAL)
    PRIMARY KEY (uniprot_id, source)
```

`source ∈ {"uniprot", "interpro", "hgnc"}`. Payload schema is stable
across versions — see `_parse_uniprot`, `_parse_interpro`,
`fetch_hgnc_to_uniprot` in `enrichment.py`.

---

## Condition vocabulary

Single source of truth is `CONDITIONS` in `streamlit_app.py`. Codes:

```
day6, day9, hs, q35, q40,
myosin_ts_15, myosin_ts_25,
paramyosin_ts_15, paramyosin_ts_25
```

Original Table S1 column names (for anyone rebuilding from Excel):

| Code               | Table S1 average column                                 |
|--------------------|---------------------------------------------------------|
| `day6`             | `AvgLog₂(wt day6/wt day1).conformation`                 |
| `day9`             | `AvgLog₂(wt day9/wt day1).conformation`                 |
| `hs`               | `AvgLog₂(wt heat-shock 35°C/wt 20°C).conformation`      |
| `q35`              | `AvgLog₂(Q35/Q24).conformation`                         |
| `q40`              | `AvgLog₂(Q40/Q24).conformation`                         |
| `myosin_ts_15`     | `AvgLog₂(myosin-ts 15°C/wt 15°C).conformation`          |
| `myosin_ts_25`     | `AvgLog₂(myosin-ts 25°C/wt 25°C).conformation`          |
| `paramyosin_ts_15` | `AvgLog₂(paramyosin-ts 15°C/wt 15°C).conformation`      |
| `paramyosin_ts_25` | `AvgLog₂(paramyosin-ts 25°C/wt 25°C).conformation`      |

---

## Recipes

**Load a protein's residue values for coloring a structure:**

```python
import pandas as pd
per_res = pd.read_parquet("data/lipms_per_residue_agg.parquet")
vit6 = per_res.query("uniprot_id == 'P18948' and condition == 'day6'")
res2fc = dict(zip(vit6['residue'], vit6['log2fc_max_abs']))
```

**Feature matrix for ML (predict residue-level aging change from structure):**

```python
long = pd.read_parquet("data/lipms_per_residue_long.parquet")
# Feature engineering per (residue, condition):
X = long.assign(
    is_semi   = (long.peptide_type == 'semi').astype(int),
    is_half   = (long.peptide_type == 'half').astype(int),
    is_full   = (long.peptide_type == 'full').astype(int),
    is_idr    = (long.idr_annotation == 'idr').astype(int),
    ss_helix  = (long.secondary_structure == 'H').astype(int),
    ss_sheet  = (long.secondary_structure == 'E').astype(int),
    has_domain= long.interpro_id.notna().astype(int),
)[['uniprot_id','condition','residue','avg_residue_depth',
   'is_semi','is_half','is_full','is_idr','ss_helix','ss_sheet','has_domain']]
y = long['log2fc']
```

Join `protein_features.parquet` on `uniprot_id` for protein-level covariates
(degree centrality, chaperone-client status, tissue specificity).

**Resolve worm → human, then to a UniProt for cross-species mapping:**

```python
import json, pandas as pd
from enrichment import fetch_hgnc_to_uniprot
o = pd.read_parquet("data/orthologs_and_disease.parquet")
row = o.query("worm_uniprot == 'P09446'").iloc[0]
for x in json.loads(row['orthologs']):
    info = fetch_hgnc_to_uniprot(x['hgnc_id'])
    print(x['hgnc_id'], '→', info['accession'], info['gene_name'], info['protein_name'])
```

---

## App entry points

| Entry point                             | Purpose |
|-----------------------------------------|---------|
| `streamlit_app.py`                      | Streamlit UI — protein page with structure, tracks, cards, orthologs. |
| `enrichment.py::fetch_uniprot`          | Cached UniProt REST lookup. |
| `enrichment.py::fetch_interpro_live`    | Cached InterPro REST lookup. |
| `enrichment.py::fetch_hgnc_to_uniprot`  | HGNC → reviewed human UniProt via REST search. |
| `enrichment.py::interpro_local`         | Offline fallback from `interpro_domains.parquet`. |

Streamlit caches wrap the enrichment functions with a one-week TTL
(`cached_uniprot`, `cached_interpro_live`, `cached_hgnc`). The SQLite
cache underneath is TTL-free — clear it manually to force refetch.

---

## What each downstream module gets, for free

### ML module (residue-level aging predictor)

Read `lipms_per_residue_long.parquet` + join `protein_features.parquet` on
`uniprot_id`. Target: `log2fc` (regression) or `abs(log2fc) > 1 AND
adj_pval < 0.05` (binary "significantly changed"). Features already
available:

- **Structural**: `avg_residue_depth`, `secondary_structure`, `idr_annotation`.
- **Sequence context**: `peptide_type`, `pk_cut_residue`.
- **Domain context**: `interpro_id` (join intervals from `interpro_domains.parquet`).
- **Protein context**: `mol_weight`, `charge`, `log2_empai`, `degree_centrality`,
  `hsp90_client`, `tric_client`, `tissue_specificity`.
- **Extend easily by adding**: pLDDT (parse from PDB `B-factor` field of
  the local AlphaFold files), SASA (Bio.PDB DSSP), backbone φ/ψ,
  distance to nearest cysteine, distance to closest active site
  (`domain_type == 'active_site'` in InterPro).

### Conversational chat module

- **Structured retrieval slot**: given a UniProt or gene, `orthologs_and_disease.parquet`
  + `protein_features.parquet` + `enrichment_cache.sqlite` give the LLM
  a JSON blob of everything it needs to answer "what does this protein do,
  what does it look like in aging, what human diseases is it linked to".
- **Grounded answer template**: paint the top-3 residue windows by
  `log2fc_max_abs` in the requested condition, cross-reference with
  `interpro_domains.parquet` and `idr_annotation` to say "the change is
  in the X domain" / "the change is in a disordered loop between X and Y".
- **Structure-biology reasoning slot**: for surface-vs-buried context, use
  `avg_residue_depth`; for stability context, use `secondary_structure`
  and pLDDT parsed from the PDB. These are all key-lookup, not
  hard-to-generate.

---

## Version & regeneration

**v0** (this document). Built by the plan `plan_build-interactive-lip-ms-aging-structure`
in frame `5af75820-c26c-4096-bf1c-96ce4377d911`. Every parquet is a
lineage-tracked artifact — `host.artifacts()` and
`host.lineage[<version_id>]` return the code that produced it.

To regenerate from raw:

1. Table S1 (`Peptide-level data`, header row 3) → `lipms_peptide_unified.parquet`
2. Join Table S3 Peptide features (header row 2) → adds structural columns
3. Explode peptides to residues → `lipms_per_residue_long.parquet`
4. Group by (uniprot_id, condition, residue) with max-abs pick →
   `lipms_per_residue_agg.parquet`
5. Table S3 Protein features (header row 2) → `protein_features.parquet`
6. `c_elegans_interpro_all_fields.csv` → `interpro_domains.parquet`
7. `simplemine_results_all_celegans.csv` → parse HGNC ortholog + disease
   → `orthologs_and_disease.parquet`
8. `ls UP000001940_6239_CAEEL_v2/*.pdb` → `alphafold_index.parquet`

---

## Known limitations

- **Coverage is peptide-level, not residue-level.** Median max coverage is
  11 % of a protein's residues — LiP-MS reports where it can cut, not
  everywhere the protein exists. Grey residues on the 3D structure mean
  "not observed", not "no change".
- **The 361 missing AlphaFold structures** are UniProt entries newer than
  the v2 bundle. On-demand fetch from AlphaFoldDB is straightforward if
  needed.
- **Ortholog TF caveats**: for transcription factors and other
  fast-evolving families, SimpleMine's ortholog inference can be
  misleading (e.g. SKN-1 → TCP1 seen in the demo). Always sanity-check
  TF orthologs against the primary literature.
- **Disease terms are `By Orthology`**: they come from the human
  ortholog's disease record, not from any worm phenotype. Interpret
  accordingly.
