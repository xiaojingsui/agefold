"""
public_content.py — plain-language content for the public "Explore aging" track.

Pure data + text only (no Streamlit imports) so it can be edited, versioned, and
fact-checked independently of the rendering code. Every factual claim in the
featured stories is grounded in the same retrieval context the researcher chat
uses (see app/PUBLIC_README.md for the accuracy protocol).
"""
from __future__ import annotations

# Featured proteins — curated for compelling, accurate aging + disease stories.
# Every number and disease link below is grounded in the LiP-MS dataset and the
# retrieval context (verified in Step 2). "signal" fields describe real measured
# aging changes; "human_health" is only present where a real worm→human ortholog
# with pathogenic ClinVar variants exists.
FEATURED: list[dict] = [
    {
        "uid": "P34697",
        "gene": "sod-1",
        "icon": "🛡️",
        "title": "sod-1 — the cell's rust-protection",
        "hook": "Your cells make their own antioxidant. This is the worm's version — and it's the same gene that, when broken in humans, causes ALS.",
        "story": (
            "Every time a cell burns fuel it makes 'sparks' — reactive oxygen molecules that "
            "damage everything they touch, a bit like rust. **sod-1** is the enzyme that mops "
            "up those sparks before they do harm. In our data, a stretch of this protein "
            "(around positions 82–137) measurably changes shape as the worm ages, right where "
            "the enzyme does its chemistry."
        ),
        "why_matters": (
            "Oxidative damage is one of the oldest theories of why we age. Watching the very "
            "enzyme that fights it lose its shape over a lifetime puts a molecular face on that idea."
        ),
        "human_health": (
            "The human version of this gene is **SOD1**. Mutations in it cause a form of "
            "**ALS (amyotrophic lateral sclerosis)** — the disease that affected Stephen Hawking. "
            "Strikingly, dozens of the known human disease mutations fall on the very stretch of "
            "the protein (around positions 82–137) that we see changing shape as the worm ages."
        ),
    },
    {
        "uid": "P10984",
        "gene": "act-2",
        "icon": "🧱",
        "title": "act-2 — the scaffolding inside every cell",
        "hook": "Actin is the internal scaffolding that gives cells their shape and lets them move. It's one of the most conserved proteins in all of life.",
        "story": (
            "Cells aren't just bags of fluid — they have an internal skeleton, and **actin** is "
            "its main building block. It's so important that its shape has barely changed across "
            "a billion years of evolution. In our data, several regions of worm actin shift with "
            "age, the strongest around positions 30–104."
        ),
        "why_matters": (
            "When the cell's scaffolding starts to warp, everything built on it — muscle, nerve "
            "connections, cell shape — is affected. It's a window into how aging reaches the "
            "structural core of a cell."
        ),
        "human_health": (
            "The human version is **ACTG1/ACTB**. Mutations cause **inherited hearing loss** and "
            "**Baraitser-Winter syndrome**, a developmental disorder. The worm and human proteins "
            "are ~98% identical, so the correspondence is unusually direct."
        ),
    },
    {
        "uid": "P18948",
        "gene": "vit-6",
        "icon": "🥚",
        "title": "vit-6 — the yolk protein that changes the most",
        "hook": "One of the most dramatically age-changed proteins in the whole dataset.",
        "story": (
            "**vit-6** makes yolk — the nutrient package a mother worm loads into her eggs. As "
            "worms age they keep pumping out yolk even when they can no longer reproduce, and it "
            "piles up in the body. Our data shows enormous shape changes across this protein — a "
            "region near position 454–500 changes more than 20-fold, putting it among the very "
            "largest age-related changes we measured."
        ),
        "why_matters": (
            "This runaway yolk production is a classic example of a program that was useful in "
            "youth becoming harmful later — a leading idea for why aging happens at all."
        ),
    },
    {
        "uid": "P02566",
        "gene": "unc-54",
        "icon": "💪",
        "title": "unc-54 — the motor that powers muscle",
        "hook": "Muscle contracts because of a molecular motor called myosin. As the worm ages, its motor changes shape — and the worm slows down.",
        "story": (
            "**unc-54** is the main **myosin** in worm muscle — the protein that physically pulls "
            "muscle fibres together so the animal can move. Our data shows shape changes in "
            "several places along this long protein, including its tail region (around 1315–1380)."
        ),
        "why_matters": (
            "Loss of muscle and strength (sarcopenia) is one of the most universal features of "
            "aging across animals. Seeing the muscle motor itself remodel with age connects that "
            "everyday experience to a specific molecule."
        ),
    },
    {
        "uid": "G5EE01",
        "gene": "daf-18",
        "icon": "🚦",
        "title": "daf-18 — a brake on growth and a guard against cancer",
        "hook": "This gene tells cells when to stop growing. In humans, losing it is one of the most common events in cancer.",
        "story": (
            "**daf-18** is a molecular brake: it counteracts the growth signals that tell cells to "
            "divide and take up nutrients. In worms it's also part of the network that controls "
            "lifespan. Our data shows several regions changing shape with age, including one near "
            "position 841–854."
        ),
        "why_matters": (
            "The same pathway that sets a worm's lifespan sets a cell's growth limits. This is one "
            "of the clearest bridges between the biology of aging and the biology of cancer."
        ),
        "human_health": (
            "The human version is **PTEN**, one of the most important tumour-suppressor genes. "
            "Inherited mutations cause **PTEN hamartoma tumour syndrome** and **Cowden syndrome**, "
            "which raise cancer risk. Over a thousand pathogenic human variants map onto this worm "
            "protein."
        ),
    },
    {
        "uid": "P54812",
        "gene": "cdc-48.2",
        "icon": "♻️",
        "title": "cdc-48.2 — the cell's recycling machine",
        "hook": "This protein pulls damaged proteins apart so they can be recycled. When it fails in humans, it causes dementia and ALS.",
        "story": (
            "Cells constantly need to clear out broken or misfolded proteins. **cdc-48.2** is a "
            "powerful molecular machine that grabs damaged proteins, unfolds them, and hands them "
            "to the cell's disposal system. Our data shows shape changes in several regions as the "
            "worm ages, including near position 262–270."
        ),
        "why_matters": (
            "A decline in protein quality-control — the cell's ability to clear its own garbage — "
            "is a hallmark of aging. This is the machine at the heart of that clean-up."
        ),
        "human_health": (
            "The human version is **VCP/p97**. Mutations cause a devastating combination of "
            "**frontotemporal dementia and ALS**, as well as a form of **Charcot-Marie-Tooth** "
            "nerve disease. The disease mutations cluster in this protein's working parts."
        ),
    },
    {
        "uid": "P09446",
        "gene": "hsp-1",
        "icon": "🔧",
        "title": "hsp-1 — the cell's repair crew",
        "hook": "This is a chaperone — a protein whose whole job is to help other proteins fold correctly and rescue them when they go wrong.",
        "story": (
            "Proteins have to fold into precise 3D shapes to work. **hsp-1** is a **chaperone**: it "
            "holds other proteins while they fold, refolds ones that have come undone, and keeps "
            "the whole system from gumming up. Our data shows its own shape shifting with age in "
            "four regions, including near position 541–551."
        ),
        "why_matters": (
            "As we age, the repair crew itself starts to falter — and when the folding helpers "
            "lose their edge, damaged proteins accumulate everywhere. Watching the chaperone age "
            "is watching the cell's maintenance system wear down."
        ),
    },
]

# Plain-language glossary shown in the "Aging science basics" expander.
GLOSSARY: list[tuple[str, str]] = [
    ("Protein",
     "A tiny molecular machine. Your body is built and run by hundreds of thousands of different proteins — they digest food, carry oxygen, contract muscle, and read your DNA."),
    ("Protein shape (structure)",
     "A protein only works if it folds into the right 3D shape, like a key cut to fit a lock. Change the shape and you change — or break — what it does."),
    ("C. elegans",
     "A tiny (1 mm) roundworm used in labs worldwide. It ages through its whole life in about two to three weeks, which lets scientists watch aging from start to finish quickly."),
    ("Aging (in this project)",
     "We compare young adult worms to older worms and ask: which proteins have changed their shape? Those shape changes are the fingerprints of aging at the molecular level."),
    ("AlphaFold",
     "An AI system that predicts a protein's 3D shape from its genetic sequence. The blobs and ribbons you can spin around in this app are AlphaFold's predicted structures."),
    ("Ortholog",
     "The 'same' gene in a different species — for example, a worm gene and its human cousin that descend from a shared ancestor and do a similar job."),
    ("Pathogenic variant / mutation",
     "A change in a gene's spelling that is known to cause disease in people. We look at where human disease mutations fall on the worm protein's shape."),
    ("Chaperone",
     "A helper protein whose job is to fold other proteins correctly and rescue them when they go wrong — a kind of cellular repair crew."),
    ("Antioxidant",
     "A molecule that neutralises the reactive 'sparks' (free radicals) made when cells burn fuel, before they can damage the cell."),
    ("Conformational change",
     "A fancy way of saying 'a change in shape.' It's the central thing this dataset measures for thousands of proteins as the worm ages."),
]

# Landing-page copy.
# Quick-start chips for the chat-forward public landing: (uid, gene, blurb).
# Curated to proteins with rich featured stories so the first chat lands well.
QUICKSTART = [
    ("P34697", "sod-1",  "An antioxidant enzyme linked to human ALS"),
    ("P09446", "hsp-1",  "A chaperone that keeps other proteins in shape"),
    ("P18948", "vit-6",  "A yolk protein that changes dramatically with age"),
    ("G5EE01", "daf-18", "A longevity gene — the worm's version of human PTEN"),
]

LANDING = {
    "headline": "See how aging reshapes the machines of life",
    "intro": (
        "Your body is run by tiny molecular machines called **proteins**. Each one only works "
        "if it folds into just the right shape. As we age, many proteins slowly lose their shape "
        "— and that may be one of the deep reasons our bodies change over time.\n\n"
        "This project measured the shape of **thousands of proteins** across the life of a tiny "
        "worm, *C. elegans*, which ages completely in just a few weeks. Below, explore some of "
        "the most striking stories — and spin the real 3D structures yourself."
    ),
    "why_worm": (
        "**Why a worm?** *C. elegans* is one of biology's most powerful tools for studying aging. "
        "It's transparent, only about a millimetre long, and lives its entire life in two to "
        "three weeks — so researchers can watch aging unfold from birth to old age in a fortnight. "
        "Remarkably, most of its genes have close human cousins, so what we learn in the worm "
        "often points directly to human biology."
    ),
    "dataset_line": "6,823 proteins tracked across the worm's whole life.",
}
