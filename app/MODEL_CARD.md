# Model Card — LiP-MS Aging Vulnerability Predictor

**Current version: v2 (baseline features + ESM-2 35M PCA-32).** v1 (baseline only) is kept as fallback and documented for comparison.

**Model**: LightGBM binary classifier + isotonic calibration
**Target**: `y = 1` if a residue has |log2FC| > 1 AND adj-p < 0.05 in **day 6 OR day 9** aging conditions (relative to day 1); else `y = 0`.
**Prediction unit**: one probability per (protein, residue).
**Training frame**: 5,252 C. elegans proteins × their observed residues = 382,980 examples. Positive class rate 12.7 % (48,704 positive residues).

## v2 changes

Two feature groups added:
1. **Structural features** (SASA, Cα contact number ≤8Å, distance to InterPro active/binding sites) — computed with Biopython + custom distance maths on the 5,015 local AlphaFold PDBs.
2. **ESM-2 35M embeddings** (per-residue, layer 12, 480-dim) reduced to top 32 PCs via IncrementalPCA (48% variance retained). Sequences >1022 aa chunked with 128-aa overlap; overlapping positions averaged.

**Ablation result** (5-fold GroupKFold by protein):

| Model | n features | OOF AUC | OOF AP | ΔAUC vs baseline |
|-------|-----------|---------|--------|-------------------|
| baseline (v1) | 28 | 0.6789 | 0.3018 | +0.0000 |
| + structural | 32 | 0.6625 | 0.2740 | **−0.0164** |
| + ESM | 60 | **0.6994** | **0.3308** | **+0.0205** |
| all (baseline + structural + ESM) | 64 | 0.6966 | 0.3294 | +0.0177 |

**v2 is `baseline + ESM`.** Adding structural features *hurt* the model — SASA and contact number carry no signal beyond what pLDDT + avg_residue_depth (already in the baseline) provided, and their per-fold variance made the fit worse. This is a robust finding across all 5 folds. Structural features are computed and saved but not used in v2.

**Recovery of v1's known false negatives** (proteins v1 completely missed at p>0.5):

| Protein | Observed positives | v1 p>0.5 | v2 p>0.5 | v1 p>0.3 | v2 p>0.3 |
|---------|--------------------|---------:|---------:|---------:|---------:|
| HSP-1 (P09446) | 82 | 0 | 2 | 27 | **86** |
| ACT-2 (P10984) | 57 | 0 | 1 | 9 | **96** |

v2 still under-predicts high-precision hits on these low-abundance proteins because global isotonic calibration is dominated by the majority (high log2_empai) proteins. But at the relative-ranking level (p>0.3), v2 recovers substantial signal that v1 missed.

**v2 top-percentile precision**:

| Score threshold | Residues (of 382,980) | Observed positive rate |
|----|-----:|-----:|
| top 1 % (score ≥ 0.940) | 3,831 | **90.7 %** (v1: 88.6 %) |
| top 5 % (score ≥ 0.374) | 20,937 | **64.2 %** (v1: 52.8 %) |
| top 10 % (score ≥ 0.251) | 38,728 | 52.8 % (v1: 45.8 %) |

---

## Intended use

Prioritize residues in a protein for follow-up work — e.g. which region of a hit protein to inspect first on the 3D structure, or which position to mutate for a stability assay. The score is meant as a **ranking**, not a definitive label.

## Not intended use

Do **not** treat a high score as a claim that a residue plays a functional role in aging. Do **not** treat a low score as evidence the residue is stable — the model has known coverage gaps (see failure modes).

---

## Performance (5-fold GroupKFold by protein — v2)

| Metric | v1 (baseline) | **v2 (baseline + ESM)** |
|--------|--------------:|--------:|
| Out-of-fold AUC-ROC | 0.679 | **0.699** |
| Out-of-fold Average Precision | 0.302 | **0.331** |
| Fold-to-fold AUC range | 0.643–0.712 | 0.668–0.732 |
| Fold-to-fold AUC σ | 0.024 | 0.023 |
| Calibration after isotonic | ≤0.02 MAE | ≤0.02 MAE |
| Top-1 % score → observed positive rate | 88.6 % | **90.7 %** |
| Top-5 % score → observed positive rate | 52.8 % | **64.2 %** |

Baseline positive class rate is 0.127; v2 AP is a 2.6× lift.

Splits are grouped by `uniprot_id` — no protein appears in both train and test. Numbers reflect what to expect on a **held-out, previously unseen protein**.

See `cv_diagnostics.png` for ROC, PR, calibration curves and top-15 feature importances.

## Features (28)

**Sequence context** (4): peptide type (semi/half/full one-hot), PK cut residue hydrophobic (F/L/I/V/W/Y).

**Structural** (7): average residue depth, pLDDT (from local AlphaFold), IDR flag, secondary-structure one-hot (H/E/T/-).

**Domain** (3): has_domain (InterPro), in_active_site, in_binding_site.

**Position** (2): relative position within protein, number of peptides covering this residue.

**Protein-level** (12): seq_length, molecular weight, charge, log2(emPAI), degree centrality, HSP90 client, TRiC client, tissue-specific, pro-longevity, anti-longevity, age-dependent ribosome pausing, ubiquitinated.

Missing values are handled natively by LightGBM (NaN goes to a learned split direction).

Top feature by gain: `n_peptides` — the model implicitly learns "how much this protein was measured" as its strongest signal.

## Known failure modes

1. **Protein-level features dominate.** Six of the top-7 features by gain are protein-level (`n_peptides`, `log2_empai`, `charge`, `seq_length`, `degree_centrality`, `mol_weight`). Residue-level structural features (rel_pos, pLDDT, avg_residue_depth) rank lower. This means the model behaves largely as a **protein ranker with residue nuance in high-signal proteins**, not a true residue-level model.
2. **False negatives on low-coverage proteins.** Preview on HSP-1 (P09446, 82 observed positive residues) and ACT-2 (P10984, 57 observed positive residues) shows **zero residues with p > 0.5** despite real signal — the low log2_empai and n_peptides for these proteins push scores below threshold. Do not conclude a protein has no vulnerable residues from a low score.
3. **True positives cluster in high-scoring proteins.** VIT-6 (P18948) has 1,120/1,129 observed positives predicted p > 0.5 — the model recognises "this is a highly-changing protein" but doesn't strongly discriminate within it. Fine-grained ranking within a single high-signal protein is not trustworthy.
4. **Underrepresented proteins.** Proteins with n_peptides < 5 or in rare tissues are underrepresented in the training set; predictions there should be treated as noise.
5. **Missing AlphaFold coverage.** 361 measured proteins have no local AF-v2 structure (TrEMBL entries newer than the UP000001940 bundle). Their pLDDT and residue depth are NaN; the model still runs but is missing one of the useful residue-level signals.
6. **Aging-specific target.** The binary label uses **day 6 or day 9 vs day 1 only**. Residues that change in heat shock or polyQ but not aging are labelled negative — even though they may be biologically interesting. Use the observed heatmap for a full-conditions view.
7. **Peptide-averaged residue values.** Where multiple peptides cover a residue, the target uses the maximum |log2FC| across them. This is a heuristic — a residue observed by conflicting peptides gets whichever peptide flagged strongest change.

## Data & code

- Training matrix: `data/training_matrix.parquet` (382,980 rows × 33 columns)
- Predictions: `data/residue_predictions.parquet` (382,980 rows)
- Model: `data/model_full.txt` (LightGBM native format, 36 trees)
- Calibrator: `data/isotonic_calibrator.pkl` (scikit-learn IsotonicRegression)
- Inference helper: `app/predictor.py` — `score_uniprot(uid)`, `top_residues(uid, k)`, `score_from_matrix(X)`
- Training / CV code: reproducible from the cells in this session's frame `5af75820-c26c-4096-bf1c-96ce4377d911`

## Roadmap

- Add per-residue structural context (SASA, distance to active site, distance to interface, contact number) — currently only depth + pLDDT are used
- Add a sequence-context model (protein language model embedding per residue) — likely narrows the protein-level dominance
- Multi-task heads for the other 7 conditions (stress / polyQ / ts-mutant) — currently only aging is targeted
