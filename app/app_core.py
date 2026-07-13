"""
app_core.py — shared foundation for the LiP-MS Aging Structure Viewer.

Single source of truth for:
  • the visual theme (Arial, teal-on-white Human Proteostasis Network house
    style) and small HTML helpers (hero, stat pills)
  • the 9-condition vocabulary
  • every cached data loader (parquets, AlphaFold index, live enrichment)
  • the diverging colormap + per-residue color function used by py3Dmol
  • build_ctx(): the resource bundle every page module reads from

The router shell (streamlit_app.py) and all page modules import from here so
loaders are cached once and the look is consistent across Home / Researcher /
Public / Help.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
from matplotlib.colors import LinearSegmentedColormap

APP_DIR = Path(__file__).parent
DATA_DIR = APP_DIR / "data"

import sys as _sys
if str(APP_DIR) not in _sys.path:
    _sys.path.insert(0, str(APP_DIR))
import enrichment  # noqa: E402

# ---- Origami crane logo (AgeFold) ----
# Faceted paper crane in the sky-blue palette; scales crisply as inline SVG.
# `_crane_svg(px)` returns the markup sized to `px` for use in brand/header marks.
_CRANE_PATHS = (
    '<polygon points="52,54 92,20 74,52" fill="#a9d4f5"/>'
    '<polygon points="30,56 60,10 74,56" fill="#7db8ec"/>'
    '<polygon points="60,10 74,56 66,54" fill="#5a92cd"/>'
    '<polygon points="30,56 74,56 46,78" fill="#7db8ec"/>'
    '<polygon points="34,58 8,34 40,56" fill="#5a92cd"/>'
    '<polygon points="8,34 2,40 16,42" fill="#7db8ec"/>'
)


def _crane_svg(px: int = 40) -> str:
    """Return the origami-crane logo as inline SVG sized to `px` pixels.

    A bold side-profile folded crane (tsuru): wing up, tail swept right, neck and
    head to the left. Solid two-tone facets so it stays legible at small sizes.
    """
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 92" '
        f'width="{px}" height="{int(px*0.92)}" style="display:block">'
        f'<g stroke="#2f5b8f" stroke-width="2.2" stroke-linejoin="round" '
        f'stroke-linecap="round">{_CRANE_PATHS}</g></svg>'
    )

# ---------------- Config ----------------
CONDITIONS = [
    ("day6",              "WT day 6",              "aging"),
    ("day9",              "WT day 9",              "aging"),
    ("hs",                "Heat shock",            "stress"),
    ("q35",               "Q35 (polyQ)",           "polyQ"),
    ("q40",               "Q40 (polyQ)",           "polyQ"),
    ("myosin_ts_15",      "myosin-ts 15°C",        "ts-mutant"),
    ("myosin_ts_25",      "myosin-ts 25°C",        "ts-mutant"),
    ("paramyosin_ts_15",  "paramyosin-ts 15°C",    "ts-mutant"),
    ("paramyosin_ts_25",  "paramyosin-ts 25°C",    "ts-mutant"),
]
COND_LABEL = {c: l for c, l, _ in CONDITIONS}
COND_GROUP = {c: g for c, _, g in CONDITIONS}

# Curated example proteins shown throughout the app.
EXAMPLE_PROTEINS = [
    ("P18948", "vit-6",  "Yolk protein — massive aging destabilization"),
    ("P34697", "sod-1",  "Superoxide dismutase — human SOD1 / ALS link"),
    ("P09446", "hsp-1",  "HSC70 chaperone — proteostasis hub"),
    ("P41988", "skn-1",  "Stress-response transcription factor"),
]


# ================= Theme =================
_THEME_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&family=IBM+Plex+Mono:wght@400;500&display=swap');
/* ============================================================================
   Minimalist scientific visualization theme.
     • Inter throughout, strict modular type scale (one ladder, no free sizes)
     • Pure-white canvas, light-gray cards, ONE muted slate-blue accent
     • Hairline 1px borders + ultra-light shadows (no heavy chrome)
     • Generous whitespace so data leads and the UI never feels cluttered
   The blue->red per-residue data colormap is intentionally NOT themed.
   ========================================================================= */
:root {
  /* Force light rendering of native form controls regardless of the viewer's
     OS dark/light setting — keeps inputs/dropdowns light on a dark-mode Mac. */
  color-scheme: light;
  /* Type scale (rem ladder — the ONLY font sizes used) */
  --fs-xs:0.75rem; --fs-sm:0.875rem; --fs-base:1rem; --fs-md:1.125rem;
  --fs-lg:1.25rem; --fs-xl:1.5rem; --fs-2xl:1.875rem; --fs-3xl:2.25rem; --fs-hero:3rem;
  --fw-normal:400; --fw-medium:500; --fw-semi:600; --fw-bold:700; --fw-black:800;
  /* Radius */
  --r-sm:0.5rem; --r-md:0.75rem; --r-lg:1rem; --r-xl:1.25rem; --r-pill:9999px;
  /* Neutral slate text */
  --ink:#1f2933; --ink2:#3e4c59; --muted:#7b8794; --muted2:#9aa5b1;
  /* Hairline borders */
  --line:#e8ebef; --line2:#dce0e6;
  /* Single muted slate-blue accent */
  --brand:#4f6d9a; --brand-dark:#3c5680; --accent:#4f6d9a; --teal:#b9c7dd;
  /* Surfaces */
  --canvas:#ffffff; --card:#ffffff; --panel:#f4f6f9; --panel2:#f7f8fa;
  --skyblue:#ffffff; --hero-ink:#1f2933;
  /* Ultra-light shadows */
  --shadow-sm:0 1px 2px rgba(17,24,39,0.04);
  --shadow-md:0 2px 10px rgba(17,24,39,0.06);
  --shadow-lg:0 8px 28px rgba(17,24,39,0.07);
}
/* ---- Base typography ---- */
html, body { color-scheme: light; }
html, body, [class*="css"], .stApp, .stMarkdown, p, span, div, label, input, textarea, select, button {
  font-family:'Inter',system-ui,-apple-system,'Segoe UI',Roboto,sans-serif !important;
  -webkit-font-smoothing:antialiased;
}
/* The Inter override above breaks Streamlit's Material Symbols icon ligatures,
   so the expander toggle renders as literal "arrow_right" text. Hide it — the
   expander header label is still fully clickable. */
[data-testid="stExpanderIcon"],
summary [data-testid="stIconMaterial"],
details summary span[data-testid="stExpanderIcon"] { display:none !important; }
/* The same Inter override makes chat-message avatars render their raw icon
   ligature ("face", "smart_toy") instead of the glyph. Hide the avatars — the
   message bubbles read fine without them. */
[data-testid="stChatMessageAvatarUser"],
[data-testid="stChatMessageAvatarAssistant"],
[data-testid="stChatMessageAvatarCustom"] { display:none !important; }
body, p, .stMarkdown, .stMarkdown p { color:var(--ink2); font-size:var(--fs-base); line-height:1.65; }
h1,h2,h3,h4,h5 { font-family:'Inter',system-ui,sans-serif !important; color:var(--ink);
  letter-spacing:-0.02em; line-height:1.25; margin:0.3rem 0; }
h1 { font-size:var(--fs-2xl); font-weight:var(--fw-black); }
h2 { font-size:var(--fs-xl);  font-weight:var(--fw-bold); }
h3 { font-size:var(--fs-lg);  font-weight:var(--fw-semi); border-bottom:1px solid var(--line);
     padding-bottom:0.5rem; margin-top:2rem; margin-bottom:0.9rem; }
h4 { font-size:var(--fs-md);  font-weight:var(--fw-semi); color:var(--ink2); }
code, kbd, .mono { font-family:'IBM Plex Mono',ui-monospace,monospace !important; font-size:0.9em;
  background:var(--panel); padding:1px 5px; border-radius:4px; }
small, .caption, [data-testid="stCaptionContainer"] { font-size:var(--fs-sm) !important; color:var(--muted) !important; }

/* ---- App chrome: hide framework menus/footer/toolbar ---- */
#MainMenu, footer, div[data-testid="stToolbar"], div[data-testid="stDecoration"],
div[data-testid="stStatusWidget"], span[data-testid="stStatusWidget"],
button[kind="header"], [data-testid="stMainMenu"] { display:none !important; }
header[data-testid="stHeader"], .stAppHeader { background:transparent !important; box-shadow:none !important;
  height:0 !important; min-height:0 !important; }
.stApp, section[data-testid="stAppViewContainer"], [data-testid="stMain"] { background:var(--canvas) !important; }

/* ---- Generous whitespace: roomy container + gaps between blocks ---- */
.block-container { padding:2.6rem 2.4rem 4rem !important; max-width:1500px; }
[data-testid="stVerticalBlock"] { gap:1.15rem; }
[data-testid="stHorizontalBlock"] { gap:1.2rem; }
[data-testid="stMain"] [data-testid="element-container"] { margin-bottom:0.1rem; }

/* ---- Top navigation (rendered on every page via app_core.top_nav) ---- */
.topnav-brand { display:flex; align-items:center; gap:10px; font-weight:var(--fw-bold);
  font-size:var(--fs-lg); color:var(--ink); letter-spacing:-0.02em; padding-top:6px; }
.topnav-brand .tn-mark { display:inline-flex; align-items:center; }
.topnav-brand .tn-tag { font-weight:var(--fw-medium); font-size:var(--fs-sm); color:var(--muted); }
div[data-testid="stPageLink"] a {
  border-radius:var(--r-sm) !important; padding:10px 18px !important; font-weight:var(--fw-semi) !important;
  font-size:var(--fs-md) !important; color:var(--ink) !important; transition:all .15s ease; justify-content:center; }
div[data-testid="stPageLink"] a:hover { background:var(--panel) !important; color:var(--brand-dark) !important; }
div[data-testid="stPageLink"] a p { color:inherit !important; font-weight:var(--fw-semi) !important;
  font-size:var(--fs-md) !important; }
.topnav-rule { border:none; border-top:1px solid var(--line); margin:10px 0 26px 0; }

/* ---- Hero (in-page header) — light, minimalist, hairline border ---- */
.app-hero {
  background:var(--panel2); border:1px solid var(--line); border-radius:var(--r-lg);
  padding:22px 28px; margin:0 0 28px 0; color:var(--ink);
  box-shadow:var(--shadow-sm); display:flex; align-items:center; gap:18px; }
.app-hero .mark { font-size:var(--fs-xl); line-height:1; background:var(--card);
  border:1px solid var(--line); border-radius:var(--r-md); padding:10px 13px; }
.app-hero .title { font-size:var(--fs-xl); font-weight:var(--fw-bold); margin:0; color:var(--ink); letter-spacing:-0.02em; }
.app-hero .subtitle { font-size:var(--fs-sm); margin:5px 0 0 0; font-weight:var(--fw-normal); color:var(--muted); }
.app-hero .badge { margin-left:auto; background:var(--card); border:1px solid var(--line2);
  border-radius:var(--r-pill); padding:6px 16px; font-size:var(--fs-xs); font-weight:var(--fw-semi);
  white-space:nowrap; color:var(--brand); letter-spacing:0.02em; }

/* ---- Big landing hero (Home) — light card, dark text, accent word ---- */
.land-hero {
  background:var(--panel2); border:1px solid var(--line); border-radius:var(--r-xl);
  padding:64px 56px; margin:6px 0 40px 0; color:var(--ink);
  box-shadow:var(--shadow-sm); position:relative; overflow:hidden; }
.land-hero .eyebrow { font-size:var(--fs-xs); letter-spacing:0.2em; text-transform:uppercase;
  color:var(--brand); opacity:0.9; margin:0 0 20px 0; font-weight:var(--fw-semi); }
.land-hero h1 { font-size:var(--fs-hero); line-height:1.1; font-weight:var(--fw-black); margin:0;
  color:var(--ink); letter-spacing:-0.03em; max-width:940px; }
.land-hero h1 .accent { color:var(--brand); }
.land-hero p.lede { font-size:var(--fs-md); line-height:1.65; margin:24px 0 0 0; max-width:680px;
  color:var(--muted); font-weight:var(--fw-normal); }

/* ---- Sidebar (home for all controls/filters) ---- */
section[data-testid="stSidebar"] { background:var(--panel2); border-right:1px solid var(--line); }
section[data-testid="stSidebar"] .block-container { padding-top:1.6rem; padding-left:1.3rem; padding-right:1.3rem; }
section[data-testid="stSidebar"] [data-testid="stVerticalBlock"] { gap:0.9rem; }
section[data-testid="stSidebar"] label { font-size:var(--fs-sm) !important; color:var(--ink2) !important; font-weight:var(--fw-medium); }
.sb-brand { display:flex; align-items:center; gap:11px; padding:2px 2px 16px 2px;
  border-bottom:1px solid var(--line); margin-bottom:18px; }
.sb-brand .sb-mark { width:36px; height:36px; border-radius:var(--r-md); flex:0 0 auto;
  background:var(--brand); display:flex; align-items:center; justify-content:center; font-size:var(--fs-lg); }
.sb-brand .sb-name { font-size:var(--fs-base); font-weight:var(--fw-bold); letter-spacing:-0.01em; line-height:1.2; }
.sb-brand .sb-tag { font-size:var(--fs-xs); color:var(--muted); }
.sb-label { font-size:var(--fs-xs); font-weight:var(--fw-semi); text-transform:uppercase;
  letter-spacing:0.08em; color:var(--muted); margin:20px 0 6px 0; }
section[data-testid="stSidebar"] div[role="radiogroup"] { gap:3px; }

/* ---- Tabs (understated segmented) ---- */
.stTabs [data-baseweb="tab-list"] { gap:6px; border-bottom:1px solid var(--line); padding-bottom:6px; }
.stTabs [data-baseweb="tab"] { height:auto; padding:9px 18px; border-radius:var(--r-sm);
  background:transparent; font-weight:var(--fw-medium); font-size:var(--fs-sm); color:var(--muted);
  border:1px solid transparent; }
.stTabs [data-baseweb="tab"]:hover { background:var(--panel); color:var(--ink2); }
.stTabs [aria-selected="true"] { background:var(--panel); color:var(--brand-dark) !important; border:1px solid var(--line); }

/* ---- Buttons ---- */
.stButton>button { border-radius:var(--r-sm); border:1px solid var(--line2); background:var(--card);
  color:var(--ink2); font-weight:var(--fw-semi); font-size:var(--fs-lg); padding:16px 24px;
  transition:all .15s ease; text-align:left; line-height:1.4; box-shadow:none; }
/* Streamlit renders the label inside a nested markdown container whose own
   emotion-generated font-size outranks a rule on the <button>. Target the
   button by its stable testid and force EVERY inner text node to the size. */
button[data-testid^="stBaseButton"],
button[data-testid^="stBaseButton"] *,
div[data-testid="stButton"] button,
div[data-testid="stButton"] button * {
  font-size:var(--fs-lg) !important; font-weight:var(--fw-bold) !important; line-height:1.4 !important;
  color:#0284c7 !important; }
.stButton>button:hover { background:var(--panel); border-color:var(--brand); color:var(--brand-dark);
  box-shadow:var(--shadow-sm); }
.stButton>button[kind="primary"] { border:1px solid var(--brand); background:var(--brand); color:#fff;
  font-weight:var(--fw-semi); text-align:center; }
.stButton>button[kind="primary"]:hover { background:var(--brand-dark); border-color:var(--brand-dark);
  color:#fff; box-shadow:var(--shadow-md); }
.stDownloadButton>button { border-radius:var(--r-sm); font-weight:var(--fw-medium); font-size:var(--fs-sm);
  border:1px solid var(--line2); background:var(--card); color:var(--ink2); }

/* ---- Protein title ---- */
.protein-title { font-weight:var(--fw-black); font-size:var(--fs-2xl); letter-spacing:-0.02em; margin:2px 0 0 0; color:var(--ink); }
.protein-title .acc { color:var(--muted); font-size:0.5em; font-weight:var(--fw-medium); font-family:'IBM Plex Mono',monospace !important; }

/* ---- Chips ---- */
.flag-row { display:flex; gap:8px; flex-wrap:wrap; margin:14px 0 2px 0; }
.chip { background:var(--panel); border:1px solid var(--line); border-radius:var(--r-pill);
  padding:5px 14px; font-size:var(--fs-xs); font-weight:var(--fw-medium); color:var(--brand-dark); }

/* ---- Dataframes + inputs: full width, hairline border, no boxing ---- */
div[data-testid="stDataFrame"] { border:1px solid var(--line); border-radius:var(--r-md); box-shadow:var(--shadow-sm); }
div[data-testid="stDataFrame"], [data-testid="stImage"], [data-testid="stImage"] img,
[data-testid="stPlotlyChart"], [data-testid="stVegaLiteChart"] { width:100% !important; }
[data-testid="stImage"] img { border-radius:var(--r-md); }
div[data-baseweb="input"]>div, div[data-baseweb="select"]>div { border-radius:var(--r-sm) !important;
  border-color:var(--line2) !important; }
div[data-baseweb="input"]>div:focus-within, div[data-baseweb="select"]>div:focus-within { border-color:var(--brand) !important; }
.stAlert { border-radius:var(--r-md); font-size:var(--fs-sm); border:1px solid var(--line); }
div[data-testid="stChatMessage"] { border-radius:var(--r-md); background:var(--panel2); border:1px solid var(--line); box-shadow:none; }

/* ---- Stat pills ---- */
.stat-row { display:flex; gap:14px; flex-wrap:wrap; justify-content:center; margin:6px 0 24px 0; }
.stat-pill { background:var(--card); border:1px solid var(--line); border-radius:var(--r-md);
  padding:14px 22px; min-width:100px; box-shadow:var(--shadow-sm); }
.stat-pill .stat-val { font-size:var(--fs-lg); font-weight:var(--fw-bold); color:var(--ink); letter-spacing:-0.02em; }
.stat-pill .stat-lab { font-size:var(--fs-xs); color:var(--muted); margin-top:4px; text-transform:uppercase; letter-spacing:0.06em; font-weight:var(--fw-medium); }

/* ---- In-page hero heading ---- */
.hero-head { margin:10px 0 22px 0; }
.hero-title { font-size:var(--fs-2xl); font-weight:var(--fw-black); letter-spacing:-0.02em; margin:0; color:var(--ink); }
.hero-sub { font-size:var(--fs-md); color:var(--muted); margin:10px 0 0 0; max-width:680px; line-height:1.6; }

/* ---- Feature card (Home) ---- */
.feat-card { border:1px solid var(--line); border-radius:var(--r-lg); padding:30px 26px 24px 26px;
  background:var(--card); height:100%; transition:all .18s ease; box-shadow:var(--shadow-sm); }
.feat-card:hover { border-color:var(--line2); box-shadow:var(--shadow-md); transform:translateY(-2px); }
.feat-card .fc-icon { font-size:var(--fs-2xl); }
.feat-card .fc-title { font-size:var(--fs-lg); font-weight:var(--fw-semi); margin:14px 0 8px 0; letter-spacing:-0.01em; color:var(--ink); }
.feat-card .fc-body { font-size:var(--fs-sm); color:var(--muted); line-height:1.65; }

/* ---- Links ---- */
a { color:var(--brand) !important; font-weight:var(--fw-medium); text-decoration:none; }
a:hover { color:var(--brand-dark) !important; text-decoration:underline; }

/* ---- Status pill (live / offline) ---- */
.status-pill { display:inline-flex; align-items:center; gap:8px; border-radius:var(--r-pill);
  padding:6px 16px; font-size:var(--fs-sm); font-weight:var(--fw-medium); margin:2px 0 14px 0; border:1px solid var(--line); }
.status-pill .dot { width:8px; height:8px; border-radius:50%; display:inline-block; }
.status-pill.live { background:var(--panel2); color:var(--brand-dark); }
.status-pill.live .dot { background:var(--brand); box-shadow:0 0 0 3px rgba(79,109,154,0.18); }
.status-pill.off { background:var(--panel2); color:var(--muted); }
.status-pill.off .dot { background:var(--muted2); }
</style>
"""


def inject_theme() -> None:
    """Inject the shared CSS theme. Call once per rerun, in the shell."""
    st.markdown(_THEME_CSS, unsafe_allow_html=True)


def top_nav(active: str = "") -> None:
    """Render the persistent top navigation bar on every page.

    Uses st.page_link against the st.Page objects the shell stashed in
    st.session_state['_nav_pages'], so links work from anywhere (this is what
    lets a user jump back to Home from the Researcher page). `active` is the key
    of the current page ('home'/'researcher'/'public'/'help') — reserved for
    highlighting; st.page_link auto-marks the current page.
    """
    pages = st.session_state.get("_nav_pages", {})
    # brand on the left, links on the right
    cols = st.columns([3.0, 1, 1, 1, 1, 1])
    with cols[0]:
        st.markdown(
            f"<div class='topnav-brand'><span class='tn-mark'>{_crane_svg(42)}</span>"
            "<span>AgeFold</span></div>",
            unsafe_allow_html=True)
    order = [("home", "Home"), ("researcher", "Researcher"),
             ("public", "Public"), ("discover", "Discover"), ("help", "Help")]
    for (key, label), col in zip(order, cols[1:]):
        with col:
            target = pages.get(key)
            if target is not None:
                st.page_link(target, label=label, use_container_width=True)
            else:
                st.markdown(f"<div style='text-align:center;color:#94a3b8;padding:8px'>{label}</div>",
                            unsafe_allow_html=True)
    st.markdown("<hr class='topnav-rule'>", unsafe_allow_html=True)


def app_header(subtitle: str, badge: str = "", title: str = "AgeFold") -> None:
    badge_html = f'<div class="badge">{badge}</div>' if badge else ""
    st.markdown(
        f"""
<div class="app-hero">
  <div class="mark">{_crane_svg(46)}</div>
  <div>
    <p class="title">{title}</p>
    <p class="subtitle">{subtitle}</p>
  </div>
  {badge_html}
</div>""",
        unsafe_allow_html=True,
    )


def stat_pills(pills: list[tuple[str, str]]) -> None:
    """ToolUniverse-style stat row: [(value, label), ...]."""
    cells = "".join(
        f'<div class="stat-pill"><div class="stat-val">{v}</div>'
        f'<div class="stat-lab">{l}</div></div>' for v, l in pills)
    st.markdown(f'<div class="stat-row">{cells}</div>', unsafe_allow_html=True)


# ================= Cached loaders =================
@st.cache_data(show_spinner=False)
def load_per_residue() -> pd.DataFrame:
    return pd.read_parquet(DATA_DIR / "lipms_per_residue_agg.parquet")


@st.cache_data(show_spinner=False)
def load_predictions() -> tuple[pd.DataFrame, str]:
    for name, label in [("residue_predictions_v2.parquet", "v2 (baseline + ESM-2)"),
                         ("residue_predictions.parquet", "v1 (baseline)")]:
        p = DATA_DIR / name
        if p.exists():
            return pd.read_parquet(p), label
    return pd.DataFrame(columns=["uniprot_id", "residue", "p_destabilized"]), "none"


@st.cache_data(show_spinner=False)
def load_discovery() -> pd.DataFrame:
    p = DATA_DIR / "discovery_scores.parquet"
    if p.exists():
        return pd.read_parquet(p)
    return pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_variants() -> pd.DataFrame:
    p = DATA_DIR / "variant_density.parquet"
    if p.exists():
        return pd.read_parquet(p)
    return pd.DataFrame(columns=["worm_uniprot", "worm_residue", "n_variants_all",
                                 "n_variants_pathogenic", "top_disease"])


@st.cache_data(show_spinner=False)
def load_peptides() -> pd.DataFrame:
    return pd.read_parquet(DATA_DIR / "lipms_peptide_unified.parquet")


@st.cache_data(show_spinner=False)
def load_protein_features() -> pd.DataFrame:
    return pd.read_parquet(DATA_DIR / "protein_features.parquet").drop_duplicates("uniprot_id")


@st.cache_data(show_spinner=False)
def load_interpro() -> pd.DataFrame:
    return pd.read_parquet(DATA_DIR / "interpro_domains.parquet")


@st.cache_data(show_spinner=False)
def load_af_index() -> dict:
    idx = pd.read_parquet(DATA_DIR / "alphafold_index.parquet")
    return dict(zip(idx["uniprot_id"], idx["pdb_path"]))


@st.cache_data(show_spinner=False)
def load_orthologs() -> dict:
    import json as _json
    p = DATA_DIR / "orthologs_and_disease.parquet"
    if not p.exists():
        return {}
    df = pd.read_parquet(p)
    out = {}
    for _, r in df.iterrows():
        out[r["worm_uniprot"]] = {
            "worm_gene": r["worm_gene"],
            "wormbase_id": r["wormbase_id"],
            "orthologs": _json.loads(r["orthologs"]),
            "diseases":  _json.loads(r["diseases"]),
            "concise_description": r.get("concise_description"),
        }
    return out


@st.cache_data(show_spinner=False)
def gene_to_uniprot_map(prot_feat: pd.DataFrame) -> dict:
    """Case-insensitive gene-symbol -> uniprot_id (first on collision)."""
    m = {}
    for uid, gs in zip(prot_feat["uniprot_id"], prot_feat["gene_symbol"]):
        if isinstance(gs, str):
            m.setdefault(gs.lower(), uid)
    return m


@st.cache_data(show_spinner=False)
def read_pdb(path: str) -> str:
    with open(path) as f:
        return f.read()


# Directory for structures fetched on-demand from the EBI AlphaFold API.
# Used when the pre-indexed local (OneDrive) PDB path is absent — i.e. on any
# deployment server. Kept out of git (.gitignore) and rebuilt lazily.
AF_CACHE_DIR = DATA_DIR / "af_cache"


@st.cache_data(show_spinner=False, ttl=60 * 60 * 24 * 30)
def fetch_alphafold_pdb(uid: str) -> str | None:
    """Fetch a UniProt's AlphaFold PDB from the EBI API, caching to disk.

    Returns PDB text, or None if no model exists / fetch fails. The EBI
    prediction API is queried for the current model URL (the bundled model
    version drifts over time — e.g. v2 → v6 — so we never hardcode it)."""
    import json as _json
    import urllib.request

    cache_file = AF_CACHE_DIR / f"AF-{uid}-F1.pdb"
    if cache_file.exists():
        try:
            return cache_file.read_text()
        except Exception:
            pass
    try:
        api = f"https://alphafold.ebi.ac.uk/api/prediction/{uid}"
        req = urllib.request.Request(api, headers={"User-Agent": "python-urllib"})
        with urllib.request.urlopen(req, timeout=30) as r:
            recs = _json.load(r)
        if not recs:
            return None
        pdb_url = recs[0].get("pdbUrl")
        if not pdb_url:
            return None
        with urllib.request.urlopen(pdb_url, timeout=60) as r:
            txt = r.read().decode()
        try:
            AF_CACHE_DIR.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(txt)
        except Exception:
            pass  # read-only FS is fine; st.cache_data still holds it in memory
        return txt
    except Exception:
        return None


def get_pdb_text(uid: str) -> str | None:
    """Resolve a UniProt's AlphaFold structure to PDB text.

    Order: pre-indexed local file (fast path on the author's Mac) → on-demand
    EBI fetch (deployment path). Returns None when no structure is available."""
    local = load_af_index().get(uid)
    if local and Path(local).exists():
        return read_pdb(local)
    return fetch_alphafold_pdb(uid)


@st.cache_data(show_spinner="Fetching UniProt…", ttl=60*60*24*7)
def cached_uniprot(uid: str) -> dict:
    return enrichment.fetch_uniprot(uid)


@st.cache_data(show_spinner=False, ttl=60*60*24*7)
def cached_interpro_live(uid: str) -> dict:
    return enrichment.fetch_interpro_live(uid)


@st.cache_data(show_spinner=False, ttl=60*60*24*30)
def cached_hgnc(hgnc_id: str) -> dict:
    return enrichment.fetch_hgnc_to_uniprot(hgnc_id)


# ================= Coloring =================
def diverging_color(v: float, vmax: float) -> str:
    """Blue (neg) -> light grey (0) -> red (pos), for py3Dmol residue coloring."""
    if not np.isfinite(v):
        return "#dddddd"
    v = max(-vmax, min(vmax, v))
    t = v / vmax  # in [-1, 1]
    if t >= 0:
        r = int(230 - t*(230-183))
        g = int(230 - t*(230-45))
        b = int(230 - t*(230-45))
    else:
        r = int(230 - (-t)*(230-45))
        g = int(230 - (-t)*(230-95))
        b = int(230 - (-t)*(230-165))
    return f"#{r:02x}{g:02x}{b:02x}"


LIPMS_CMAP = LinearSegmentedColormap.from_list("lipms", ["#2d5da0", "#e6e6e6", "#b73333"], N=256)
LIPMS_CMAP.set_bad("#f5f5f5")


# ================= Shared context bundle =================
@st.cache_resource(show_spinner=False)
def _ctx_singleton() -> dict:
    """Build the resource bundle once per session. Data loaders are themselves
    cached, so this mostly wires references together."""
    prot_feat = load_protein_features()
    return {
        "prot_feat": prot_feat,
        "gene2uid": gene_to_uniprot_map(prot_feat),
        "per_res": load_per_residue(),
        "peptides": load_peptides(),
        "interpro": load_interpro(),
        "af_index": load_af_index(),
        "read_pdb": read_pdb,
        "get_pdb_text": get_pdb_text,
        "diverging_color": diverging_color,
        "load_variants": load_variants,
        "load_orthologs": load_orthologs,
        "load_predictions": load_predictions,
        "load_discovery": load_discovery,
        "cached_uniprot": cached_uniprot,
        "cached_interpro_live": cached_interpro_live,
        "cached_hgnc": cached_hgnc,
        "CONDITIONS": CONDITIONS,
        "COND_LABEL": COND_LABEL,
        "COND_GROUP": COND_GROUP,
        "LIPMS_CMAP": LIPMS_CMAP,
        "app_header": app_header,
        "stat_pills": stat_pills,
    }


def build_ctx() -> dict:
    """Return the shared resource dict used by every page module."""
    return _ctx_singleton()
