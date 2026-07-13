# LiP-MS Aging Structure Viewer

Interactive per-residue exploration of TMT-LiP-MS conformational changes in
*C. elegans* across 9 conditions, painted onto local AlphaFold structures.

## Quick start

```bash
# Activate the environment
conda activate lipms-viewer

# From this directory
streamlit run streamlit_app.py
```

The app opens at `http://localhost:8501`.

## What you get

- **Search by UniProt ID or gene symbol** (case-insensitive; try `P18948`, `vit-6`, `hsp-1`, `daf-16`).
- **9-condition selector**: WT day 6 / day 9, heat shock, Q35, Q40, myosin-ts 15°C / 25°C, paramyosin-ts 15°C / 25°C.
- **Interactive 3D viewer** (py3Dmol): AlphaFold cartoon painted by per-residue log₂FC on a zero-centered diverging colormap; unmeasured residues grey.
- **Residue signal heatmap** across all 9 conditions, with the current condition highlighted.
- **InterPro domain track** below the heatmap.
- **Protein card**: gene, length, MW, log₂(emPAI), network degree centrality, longevity flag, ribosome pausing, turnover, ubiquitination, HSP90/TRiC client status, tissue specificity.
- **Peptide table** for the current condition sorted by |log₂FC|.

## Data files (in `data/`)

| File | What |
|---|---|
| `lipms_per_residue_agg.parquet` | 3.8M rows — (uniprot_id, condition, residue) → log2fc_max_abs, log2fc_mean, adj_pval, peptide_type, n_peptides. **Paints the structure.** |
| `lipms_peptide_unified.parquet` | 531k rows — (peptide × condition) with log2fc, p-value, structural features, InterPro. **Peptide table.** |
| `protein_features.parquet` | 9.5k rows — per-protein annotations from Table S3. **Protein card.** |
| `interpro_domains.parquet` | 58k rows — domain intervals. **Domain track.** |
| `alphafold_index.parquet` | 19.7k rows — uniprot_id → local PDB path. **Structure resolver.** |

See `../DATA_SCHEMA.md` (added in step 6) for full column dictionary.

## Roadmap

- **Step 4** — UniProt + InterPro live lookups (cached), enriching the protein card with function summary and GO terms.
- **Step 5** — Worm→human ortholog panel with ClinVar/UniProt variant preview.
- **Later** — ML module predicting which residues change with age from structural features; conversational chat interface over the same data.
