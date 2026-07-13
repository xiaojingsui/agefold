"""
Prompt templates for the protein-page chat.

Design goals:
- Ground every claim in the retrieved context — refuse when data is absent.
- Use short bracketed citations [S1]-[S6] tied to the sources dict; NO free citations.
- Prefer structural-biology reasoning over generic textbook statements.
- Distinguish observed LiP-MS data (S5) from ML predictions (S6) explicitly.
"""
from __future__ import annotations
import json
from typing import Any


SYSTEM_PROMPT = """You are a structural-biology + protein-aging assistant embedded in a viewer for a C. elegans TMT-LiP-MS aging dataset (Sui et al.). You answer questions about a single protein at a time, using ONLY the structured context provided in the user message.

Rules you must follow every turn:

1. Cite your sources. Every factual claim gets a bracketed citation from the sources block: [S1] UniProt, [S2] InterPro, [S3] AlphaFold, [S4] WormBase, [S5] LiP-MS dataset (this study), [S6] ML predictor (this study). Multiple citations use [S1, S5]. NO citations to sources outside that list.

2. Ground in the context. If a fact is not in the context, say "The context does not contain that information" and suggest what the user could look up. Do NOT rely on your training knowledge for protein-specific facts.

3. Distinguish evidence types clearly:
   - "Observed in LiP-MS" = measured peptide-level change (from aging_signal or lipms sections) → [S5]
   - "ML-predicted" = model score, not a measurement → [S6]. Always note that predictions are prioritization scores, not measurements.
   - Domains and function → [S1] or [S2] as appropriate.

4. Prefer specific residue numbers, domain names, and regions over vague phrases. When aging_signal.top_regions lists ranges, quote them.

5. Structural-biology interpretation is welcome — connecting an aging-changed region to the domain it lives in, or to solvent-exposed vs buried context — but only using facts from the context. If linking to disease, use the orthologs.diseases list and note the "by-orthology" nature of the inference.

6. Be concise. 3-6 short paragraphs is the target for a substantive question. Use bullets when listing residues, regions, or orthologs.

7. If the user asks a general question the retrieval doesn't answer (e.g. "how does LiP-MS work"), give a 2-3 sentence textbook answer WITHOUT citations, then redirect to the current protein.

8. When the context includes VARIANT OVERLAP data, cite it as [S7] and distinguish it clearly from LiP-MS measurements: variants are ClinVar entries on the HUMAN ortholog, mapped onto worm residues via pairwise alignment. They are homology inferences, not observations in this dataset. Say "pathogenic variants in the human ortholog map to worm residues X-Y" — never claim the variants were observed in the worm data. When alignment identity is low (<30%) or missing, note that the mapping is low-confidence.

9. When the context includes STRUCTURAL CONTEXT, cite it as [S8]. This is computed from the AlphaFold model: solvent accessibility (buried core vs exposed surface), three-state secondary structure (helix/sheet/loop), per-residue pLDDT confidence, and CA–CA distance to UniProt functional sites. Use it to reason about STABILITY: a buried-core aging change near a functional site is more likely to affect fold stability or activity than an exposed surface-loop change. ALWAYS respect the confidence flag — if a region is marked LOW CONFIDENCE / likely disordered, say the structural interpretation is uncertain there rather than over-reading it. These are AlphaFold-predicted structural properties, not experimental structures; say so when it matters. The functional-site distance is geometric (CA–CA in the monomer model).
"""


SYSTEM_PROMPT_PUBLIC = """You are a warm, encouraging science communicator embedded in a public "Explore aging" tool. It shows how proteins change shape as a tiny worm (C. elegans) ages, measured by a technique called LiP-MS (Sui et al.). You answer questions from curious members of the public — people who are NOT scientists but are genuinely interested in the biology of aging.

Your audience has no biology background. Write for a bright, interested teenager or an adult reading a good popular-science article.

Rules you must follow every turn:

1. Ground everything in the provided context — same discipline as a researcher would use. Never invent protein-specific facts. If the context doesn't contain something, say so gently ("We don't have that measurement for this protein") rather than guessing. Your accuracy must be identical to the technical version; only the language changes.

2. No unexplained jargon. Never use "log2FC", "LiP-MS", "pLDDT", "InterPro", "adj-p", "alignment identity", "conformational", "residue" without translating. Prefer: "changes shape", "a spot on the protein", "position 82", "a region", "the human version of this gene". If you must introduce a real term, define it in the same sentence in plain words.

3. Use analogies and everyday images (a key fitting a lock, scaffolding, a repair crew, rust, recycling) — but keep them accurate. An analogy should illuminate the real biology, never replace it with something false.

4. Do NOT show bracketed source codes like [S1] or [S5] to the user — they're distracting for this audience. The information is still grounded in those sources under the hood; a "Sources" panel is shown separately.

5. Be honest about uncertainty and inference. When something is a comparison to the human version of a gene rather than a direct measurement in people, say plainly "the human cousin of this gene..." and note it's based on similarity. Never overstate a disease link or imply the worm data proves something about human patients.

6. Keep it short and inviting. 2-4 short paragraphs. End on why it connects to understanding aging, when natural. Warm, curious, never condescending.

7. If asked something the tool can't answer, give a friendly 2-3 sentence general explanation, then gently steer back to what they CAN explore here.

8. Positions and disease links, when you use them, come only from the context. You may say "the biggest changes happen around position 82" if the context says so — plain numbers are fine; it's the coded jargon you avoid.
"""


CONTEXT_TEMPLATE = """<protein_context>
Protein: {gene} ({uid}) — {name}
Length: {length} aa · MW: {mw} kDa

FUNCTION
{function}
Keywords: {keywords}
GO biological process: {go_bp}
GO molecular function: {go_mf}
Subcellular: {subcellular}

PROTEIN-LEVEL ANNOTATIONS
log2(emPAI) abundance: {empai}
Network degree centrality: {centrality}
HSP90 client: {hsp90}
TRiC client: {tric}
Tissue specificity: {tissue}
Longevity flag: {longevity}
Age-dependent ribosome pausing: {ribo}
Age-dependent turnover: {turnover}
Ubiquitinated: {ubiq}

INTERPRO DOMAINS ({n_domains} total, top {shown_domains} by size)
{domains}

LIP-MS OBSERVATIONS
Total unique peptides observed: {n_peptides_total}
Per-condition peptide counts: {per_cond}

AGING SIGNAL (day6 and day9, |log2FC|>1 AND adj-p<0.05)
Significant aging peptides: {n_sig_aging} across {n_regions} contiguous regions
Top regions by max |log2FC|:
{regions}

ML PREDICTIONS (model {model_version}, scope: THIS PROTEIN ONLY)
Residues of THIS protein scored: {n_scored}; of those, {n_hi} are high-confidence (p>0.5) on THIS protein.
(These counts are per-protein, not genome-wide.)
Top 10 predicted residues on THIS protein:
{ml_residues}

HUMAN ORTHOLOGS AND DISEASE
Concise WormBase description: {desc}
Top orthologs: {orth}
Disease associations (by orthology): {diseases}

VARIANT OVERLAP (ClinVar variants in human ortholog, mapped through pairwise alignment)
{variant_overlap}

STRUCTURAL CONTEXT (AlphaFold-derived: solvent accessibility, secondary structure, pLDDT confidence, distance to functional sites)
{structural_context}

SOURCES
{sources}
</protein_context>

USER QUESTION: {question}"""


def format_context_prompt(ctx: dict[str, Any], question: str) -> str:
    """Render a protein context dict + question into the LLM user message."""
    if "error" in ctx:
        return f"UNKNOWN PROTEIN: {ctx['error']}\n\nUSER QUESTION: {question}"

    ident = ctx["identity"]
    fn = ctx["function"]
    ann = ctx["protein_annotations"]
    lip = ctx["lipms"]
    ag = ctx["aging_signal"]
    ml = ctx["ml"]
    orth = ctx["orthologs"]
    doms = ctx["domains"]
    src = ctx["sources"]
    vo = ctx.get("variant_overlap") or {}

    doms_str = "\n".join(
        f"  · {d['name']} ({d['interpro_id']}, {d['type']}): residues {d['start']}-{d['end']}"
        for d in doms
    ) or "  (none annotated)"

    per_cond = ", ".join(f"{c}: {v['n_peptides']}" for c, v in lip.get("per_condition", {}).items()) or "(no LiP-MS data)"

    def _region_str(r):
        head = f"  · residues {r['start']}-{r['end']}, {r['n_peptides']} peptides, max |log2FC|={r['max_abs_log2fc']}"
        dom = r.get("in_domains", [])
        if dom:
            tag = "; ".join(f"in {d['name']} ({d['interpro_id']}, {d['span']}, {d['overlap_aa']} aa overlap)" for d in dom)
            head += f"  [{tag}]"
        return head
    regions_str = "\n".join(_region_str(r) for r in ag.get("top_regions", [])) or "  (no residues cross the significance threshold in aging)"

    ml_str = "\n".join(
        f"  · residue {m['residue']}: p={m['p_destabilized']} (observed={m['observed']}, "
        f"observed max |log2FC|={m['max_abs_log2fc_observed']})"
        for m in ml.get("top_residues", [])
    ) or "  (no ML predictions available for this protein)"

    orth_str = ", ".join(
        f"HGNC:{o['hgnc_id']} ({o['n_methods']} methods)"
        for o in orth.get("orthologs", [])
    ) or "(no worm→human orthologs curated)"

    dis = orth.get("diseases", [])
    dis_str = "; ".join(dis[:10]) if dis else "(no disease associations)"

    sources_str = "\n".join(f"  [{k}] {v}" for k, v in src.items())

    # Structural context block
    sc = ctx.get("structural_context") or {}
    if not sc.get("available"):
        sc_str = "  (no local AlphaFold structure available for this protein)"
    else:
        pp = sc.get("per_protein", {})
        lines = [f"  Whole protein: {pp.get('pct_buried','?')}% of residues buried, median pLDDT {pp.get('median_plddt','?')}."]
        if sc.get("regions"):
            lines.append("  Per aging-region structural context:")
            for r in sc["regions"]:
                site = ""
                if r.get("near_site"):
                    site = f", near {r.get('nearest_site_type','a site')} ({r.get('min_site_distance')} Å)"
                lines.append(
                    f"    · residues {r['start']}-{r['end']}: {r['pct_buried']}% buried, "
                    f"predominantly {r['dominant_sse']}, mean pLDDT {r['mean_plddt']}"
                    f"{' [LOW CONFIDENCE / likely disordered]' if r.get('low_confidence') else ''}{site}. "
                    f"{r['interpretation']}"
                )
        sc_str = "\n".join(lines)

    # Variant overlap block
    if not vo.get("available"):
        vo_str = "  (no ortholog-based variant mapping available for this protein)"
    elif vo.get("n_pathogenic_total", 0) == 0:
        vo_str = f"  Ortholog mapped, but 0 pathogenic ClinVar variants ({vo.get('n_variants_all_total', 0)} total variants of any significance)."
    else:
        lines = [
            f"  Total pathogenic variants mapped to worm sequence: {vo['n_pathogenic_total']} "
            f"(across {vo['n_residues_pathogenic']} residues, {vo.get('n_variants_all_total', 0)} variants of any significance)",
            f"  Top diseases across the protein: {'; '.join(vo['top_diseases']) or '(none)'}",
            "  Per aging-region overlap:",
        ]
        for r in vo.get("regions", []):
            tail = f" — top diseases in region: {'; '.join(r['top_diseases_in_region'])}" if r.get("top_diseases_in_region") else ""
            lines.append(f"    · residues {r['start']}-{r['end']}: {r['n_pathogenic_in_region']} pathogenic, "
                          f"{r['n_variants_all_in_region']} total{tail}")
        vo_str = "\n".join(lines)

    def _fmt(v, blank="(unknown)"):
        if v is None: return blank
        if isinstance(v, str) and v in ("", "—"): return blank
        if isinstance(v, list):
            if len(v) == 0: return blank
            return ", ".join(str(x) for x in v[:5])
        try:
            import math
            if isinstance(v, float) and math.isnan(v): return blank
        except Exception: pass
        return str(v)

    return CONTEXT_TEMPLATE.format(
        gene=_fmt(ident.get("gene_symbol")),
        uid=ident["uniprot_id"],
        name=_fmt(ident.get("recommended_name")),
        length=_fmt(ident.get("sequence_length")),
        mw=_fmt(ident.get("molecular_weight_kda")),
        function=_fmt(fn.get("text"), blank="(no function annotation in UniProt)"),
        keywords=_fmt(fn.get("keywords")),
        go_bp=_fmt(fn.get("go_biological_process")),
        go_mf=_fmt(fn.get("go_molecular_function")),
        subcellular=_fmt(fn.get("subcellular_location")),
        empai=_fmt(ann.get("log2_empai")),
        centrality=_fmt(ann.get("degree_centrality")),
        hsp90=_fmt(ann.get("hsp90_client")),
        tric=_fmt(ann.get("tric_client")),
        tissue=_fmt(ann.get("tissue_specificity")),
        longevity=_fmt(ann.get("longevity")),
        ribo=_fmt(ann.get("age_ribosome_pausing")),
        turnover=_fmt(ann.get("age_protein_turnover")),
        ubiq=_fmt(ann.get("ubiquitination")),
        n_domains=len(doms), shown_domains=len(doms),
        domains=doms_str,
        n_peptides_total=lip.get("n_peptides_total", 0),
        per_cond=per_cond,
        n_sig_aging=ag.get("n_significant_aging_peptides", 0),
        n_regions=len(ag.get("top_regions", [])),
        regions=regions_str,
        model_version=ml.get("model_version") or "(none)",
        n_scored=ml.get("n_residues_scored", 0),
        n_hi=ml.get("n_high_confidence", 0),
        ml_residues=ml_str,
        desc=_fmt(orth.get("concise_description"), blank="(no WormBase description)"),
        orth=orth_str,
        diseases=dis_str,
        variant_overlap=vo_str,
        structural_context=sc_str,
        sources=sources_str,
        question=question,
    )
