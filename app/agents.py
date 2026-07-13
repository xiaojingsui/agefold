"""
agents.py — specialist agent registry for the multi-agent structural-proteomics
framework.

Each specialist is a focused persona over the SAME grounded retrieval context the
single-protein chat uses (retrieval.get_protein_context → [S1]-[S8]). The
orchestrator (agent_orchestrator.py) routes a question to the 1-3 relevant
specialists, runs them in parallel, and synthesizes one cited answer.

This module is pure data + light routing logic — no Streamlit, no network — so it
can be unit-tested and versioned independently.
"""
from __future__ import annotations
import re
from typing import Any

# ---------------------------------------------------------------------------
# Shared grounding rules — every specialist inherits these, then adds its lane.
# ---------------------------------------------------------------------------
_SHARED_RULES = """You are one specialist on a panel analyzing a C. elegans TMT-LiP-MS aging dataset (Sui et al.). A coordinator will merge your answer with other specialists', so STAY IN YOUR LANE and be rigorous.

Grounding rules (all specialists, every turn):
- Use ONLY the structured context in the user message. If a fact is not there, say "the context does not contain that" — never fill from training knowledge for protein-specific claims.
- Cite every factual claim with a bracketed source code from the sources block. NO citations outside the provided list.
- Distinguish evidence types: measured LiP-MS observations [S5] vs ML predictions [S6] (prioritization scores, not measurements) vs homology inferences (variants [S7], orthology-based disease) vs AlphaFold-predicted structure [S8] (not experimental).
- Prefer specific residue numbers, regions, and domain names over vague phrases.
- Be concise: 2-4 short paragraphs. You are ONE voice on a panel, not the whole answer — don't repeat context the other specialists cover; focus on your specialty.
"""

# ---------------------------------------------------------------------------
# Specialist roster
# ---------------------------------------------------------------------------
# Each agent: id, name, icon, role (one line), focus (context sections it centers
# on), sources (its citation codes), and system_prompt (shared + lane).
AGENTS: dict[str, dict[str, Any]] = {
    "structural_biologist": {
        "name": "Structural Biologist",
        "icon": "🧬",
        "role": "Fold stability, secondary structure, solvent accessibility, and proximity of aging changes to functional sites.",
        "focus": ["structural_context", "domains", "aging_signal", "function"],
        "sources": ["S8", "S2", "S3", "S1"],
        "system_prompt": _SHARED_RULES + """
YOUR LANE — Structural biology and conformational stability.
Interpret where each aging-changed region sits in the fold using the STRUCTURAL CONTEXT block [S8]: buried core vs exposed surface, secondary structure (helix/sheet/loop), pLDDT confidence, and CA–CA distance to functional sites. Connect regions to InterPro domains [S2] they fall in.
Reason about STABILITY: a buried-core change near a functional site is more likely to affect fold stability or activity than an exposed surface loop. ALWAYS respect the confidence flag — if a region is LOW CONFIDENCE / likely disordered, say the structural read is uncertain there. These are AlphaFold-predicted properties, not experimental structures; SASA/site distances are monomer geometry. Do NOT opine on disease genetics or cross-condition statistics — other specialists cover those.
""",
    },
    "aging_stress_analyst": {
        "name": "Aging & Stress Analyst",
        "icon": "📈",
        "role": "The 9-condition LiP-MS comparison — is a change aging-specific or shared with heat-shock, polyQ, or ts-mutant stress?",
        "focus": ["aging_signal", "lipms", "protein_annotations"],
        "sources": ["S5"],
        "system_prompt": _SHARED_RULES + """
YOUR LANE — Aging vs stress condition comparison (measured LiP-MS [S5]).
The dataset spans 9 conditions in 4 groups: aging (WT day6, day9), stress (heat shock), polyQ (Q35, Q40), and ts-mutant (myosin/paramyosin ts at 15/25°C). Your job: characterize the protein's conformational response ACROSS conditions. Is the aging signal (day6/day9) also seen under stress/polyQ/ts (a generic destabilization), or is it aging-SPECIFIC? Quote per-condition peptide counts and the aging top-regions with their max |log2FC|.
When DATASET ANALYSIS results are provided (cross-condition tables, aging-specific hit lists), use them and cite [S5]. Frame everything as measured observations, never predictions. Do NOT interpret 3D structure or disease genetics — other specialists handle those.
""",
    },
    "disease_genetics": {
        "name": "Disease & Genetics",
        "icon": "🩺",
        "role": "Human orthologs and disease-variant mapping — what the worm protein's human counterpart tells us, with homology caveats.",
        "focus": ["orthologs", "variant_overlap", "identity"],
        "sources": ["S4", "S7", "S1"],
        "system_prompt": _SHARED_RULES + """
YOUR LANE — Human orthologs and disease genetics.
Use the HUMAN ORTHOLOGS block [S4] and the VARIANT OVERLAP block [S7]. Variants are ClinVar entries on the HUMAN ortholog, mapped onto worm residues via pairwise alignment — they are HOMOLOGY INFERENCES, not observations in the worm data. Always say "pathogenic variants in the human ortholog map to worm residues X-Y", never "the worm has these variants". When alignment identity is low (<30%) or absent, flag the mapping as low-confidence. Note whether pathogenic variants fall inside the aging-changed regions (that overlap is the interesting cross-link), but attribute the aging change itself to [S5], not to yourself. Do NOT interpret fold stability or run cross-condition stats.
""",
    },
    "data_analyst": {
        "name": "Data Analyst",
        "icon": "📊",
        "role": "ML destabilization predictions and dataset-level / cross-protein questions from precomputed tables.",
        "focus": ["ml", "aging_signal"],
        "sources": ["S6", "S5"],
        "system_prompt": _SHARED_RULES + """
YOUR LANE — ML predictions [S6] and dataset-level analysis.
Two jobs: (1) Interpret the ML DESTABILIZATION predictions for this protein — always state these are prioritization SCORES, not measurements [S6], and note calibration/coverage caveats. (2) When DATASET ANALYSIS results are provided (cross-protein comparisons, top-movers, aging-specific hit lists, structural enrichment across all proteins), summarize them faithfully and cite [S5]/[S6] as appropriate. Give concrete numbers and rankings. Do NOT over-interpret single-residue structure or disease genetics — other specialists own those.
""",
    },
}

SPECIALIST_ORDER = ["structural_biologist", "aging_stress_analyst", "disease_genetics", "data_analyst"]

# ---------------------------------------------------------------------------
# Router — keyword signals per specialist; falls back to a default panel.
# ---------------------------------------------------------------------------
_ROUTER_KEYWORDS: dict[str, list[str]] = {
    "structural_biologist": [
        "structur", "stability", "stabil", "destabil", "fold", "buried", "exposed",
        "surface", "secondary structure", "helix", "sheet", "strand", "loop", "domain",
        "conformation", "plddt", "solvent", "accessib", "active site", "binding site",
        "disulfide", "3d", "alphafold", "sasa", "rigid", "flexib",
    ],
    "aging_stress_analyst": [
        "aging", "age", "day6", "day9", "day 6", "day 9", "stress", "heat", "heat shock",
        "polyq", "q35", "q40", "ts", "myosin", "paramyosin", "condition", "compare",
        "specific", "shared", "generic", "versus", " vs ", "temperature", "which condition",
    ],
    "disease_genetics": [
        "disease", "variant", "mutation", "pathogenic", "clinvar", "ortholog", "human",
        "als", "cancer", "genetic", "hgnc", "syndrome", "hereditary", "patient", "clinical",
        "homolog",
    ],
    "data_analyst": [
        "predict", "ml", "model", "score", "destabiliz", "confidence", "rank", "top ",
        "most", "how many", "count", "across all", "dataset", "proteome", "genome-wide",
        "list", "which proteins", "enrich", "statistic", "distribution", "compare proteins",
    ],
}

# Questions that are clearly about the whole dataset, not one protein.
_DATASET_SIGNALS = [
    "across all", "which proteins", "how many proteins", "top ", "most ", "rank",
    "dataset", "proteome", "genome-wide", "list the", "compare proteins", "overall",
    "in general", "enrich",
]


def is_dataset_level(question: str) -> bool:
    q = question.lower()
    return any(s in q for s in _DATASET_SIGNALS)


def route(question: str, max_specialists: int = 3) -> list[str]:
    """Pick the relevant specialists for a question by keyword scoring.

    Returns an ordered list of agent ids (1-max_specialists). Always returns at
    least one. Deterministic; the orchestrator may override with an LLM classify
    when scores are flat.
    """
    q = " " + question.lower() + " "
    scores: dict[str, int] = {}
    for aid, kws in _ROUTER_KEYWORDS.items():
        scores[aid] = sum(1 for kw in kws if kw in q)

    # Dataset-level questions always include the Data Analyst.
    if is_dataset_level(question):
        scores["data_analyst"] += 3

    ranked = sorted(SPECIALIST_ORDER, key=lambda a: (-scores[a], SPECIALIST_ORDER.index(a)))
    chosen = [a for a in ranked if scores[a] > 0][:max_specialists]

    if not chosen:
        # No keyword hit → sensible default panel for an open-ended protein question.
        chosen = ["structural_biologist", "aging_stress_analyst"]
    return chosen


def route_explain(question: str, max_specialists: int = 3) -> dict[str, Any]:
    """Router decision + per-specialist keyword scores, for tracing/UI."""
    q = " " + question.lower() + " "
    scores = {aid: sum(1 for kw in kws if kw in q) for aid, kws in _ROUTER_KEYWORDS.items()}
    if is_dataset_level(question):
        scores["data_analyst"] += 3
    return {"chosen": route(question, max_specialists), "scores": scores,
            "dataset_level": is_dataset_level(question)}
