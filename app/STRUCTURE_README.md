# Structural context module

Adds AlphaFold-derived **structural interpretation** to every aging hotspot:
is the changing region a buried core or an exposed surface? A helix, sheet, or
loop? Well-folded or low-confidence/disordered? Close to a functional site? This
turns *"the structure changes here"* into *"this buried, functionally coupled
region changes — the kind of change most likely to affect fold stability."*

This is source **[S8]** in the chat and the **Structural context** card in the
researcher panel.

---

## What is computed (per residue, all 6,462 measured proteins with a local AlphaFold model)

`app/data/structural_features.parquet` — 3,261,471 residues:

| column | meaning |
|--------|---------|
| `sasa` | absolute solvent-accessible surface area (Å²), biotite Shrake–Rupley (`vdw_radii="Single"`) |
| `rsa` | relative solvent accessibility = SASA / max-ASA(residue type), capped at 1.5 |
| `buried` | `rsa < 0.20` → buried core; else exposed surface |
| `sse` | three-state secondary structure from biotite `annotate_sse`: `H` helix, `E` sheet/strand, `C` loop/coil |
| `plddt` | AlphaFold per-residue confidence = CA-atom B-factor (0–100) |
| `disordered` | `plddt < 50` → very low confidence / likely flexible or disordered |

**Max-ASA normalization** uses the Tien et al. (2013) *theoretical* maximum
accessible surface areas per amino-acid type.

## Functional-site proximity

`app/data/uniprot_sites.parquet` — functional-site features fetched from UniProt
(`ft_act_site`, `ft_binding`, `ft_site`, `ft_disulfid`) for measured proteins:
**1,409 proteins with sites, 8,117 site rows** (5,444 binding · 1,522 disulfide ·
891 active · 260 other).

`app/data/site_proximity.parquet` — for each residue of a protein that has any
annotated site, the **CA–CA distance** (in the AlphaFold monomer model) to the
nearest site residue, its type, and a `near_site` flag (`≤ 8 Å`). Covers 1,341
proteins / 784,135 residues; 11.8% of residues are near a site.

## Region summary + interpretation

`structure_features.get_structural_context(uid, regions=[(start,end),...])`
returns, for each aging region: `pct_buried`, `dominant_sse`, `mean_plddt`,
`low_confidence` flag, `near_site` + `nearest_site_type` + `min_site_distance`,
and a one-line plain-language `interpretation`. The interpretation gates on
confidence first (disordered → "harder to interpret structurally"), then reads
buried-vs-surface × near-site into a stability statement.

---

## Thresholds

| threshold | value | rationale |
|-----------|-------|-----------|
| buried | RSA < 0.20 | standard buried/exposed cutoff |
| disordered / low-confidence | pLDDT < 50 | AlphaFold's own "very low" band |
| high confidence (panel dot) | pLDDT ≥ 70 | AlphaFold "confident" band |
| near functional site | CA–CA ≤ 8 Å | first coordination shell + a residue's reach |

---

## Validation

- **Numbering** aligns 1:1 with the LiP-MS data and the local AlphaFold PDBs —
  no offset (confirmed on sod-1/vit-6/act-2: every LiP-MS residue is covered).
- **SOD-1 sanity check**: all 9 UniProt-annotated Cu/Zn-binding + disulfide
  residues (68, 70, 85, 93, 102, 105, 142 + Cys79/171) come out **buried**
  (RSA < 0.12, pLDDT ~99) — the metal site really is in the β-barrel core.

## What the whole-dataset comparison shows (honest)

Across 45,671 hotspot residues vs all measured residues, aging hotspots are
**structurally similar to the background overall** (38.4% vs 38.8% buried, same
median pLDDT). Two real signals stand out:
- **1.31× enrichment near functional sites** (18.2% of hotspot residues vs 13.9%
  of background lie within 8 Å of an annotated site).
- Mild enrichment in β-sheet (18.9% vs 15.2%), mild depletion in helix.

Aging conformational change is not concentrated in exposed loops — it reaches
into ordered, functionally coupled regions somewhat more than chance.

---

## Limitations

- **AlphaFold monomer models.** SASA is computed on the single chain, so residues
  buried at a real oligomeric interface will read as "exposed"; site distances are
  monomer geometry, not the assembled complex.
- **Predicted, not experimental.** pLDDT is a confidence score, not a B-factor
  from crystallography; disorder is inferred, not measured.
- **Site coverage is sparse.** Only ~1,400 measured proteins carry any UniProt
  functional-site annotation, so `near_site` is NaN/False for most proteins — a
  bonus signal where present, not a universal feature.
- **RSA normalization is theoretical** (Tien max-ASA), not context-dependent.
- Secondary structure is AlphaFold-model-derived (`annotate_sse`), not DSSP on an
  experimental structure.

## Methods parameters

biotite 1.6.0 · `struc.sasa(..., vdw_radii="Single")` · `struc.annotate_sse` ·
Tien 2013 max-ASA · buried RSA<0.20 · disorder pLDDT<50 · near-site CA–CA ≤8 Å.
Full pass: 6,462 structures, ThreadPoolExecutor(6), ~12 min.
