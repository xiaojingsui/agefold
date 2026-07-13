"""
LiP-MS Aging Structure Viewer — top-navigation shell.

A conversational-AI resource for the C. elegans TMT-LiP-MS aging dataset. This
file is a thin router: it sets page config, injects the shared theme, and mounts
a native top navigation bar (Home / Researcher / Discover / Public / Help) via
st.navigation. All page content lives in page modules; all shared data loaders,
theme, and helpers live in app_core.

  streamlit_app.py  ── shell (this file): config + theme + top nav
  app_core.py       ── theme, config, cached loaders, colors, build_ctx()
  page_home.py      ── Home showcase / landing
  page_researcher.py→ researcher_view.render()  (full technical viewer + tabs)
  page_public.py    → public_view.render_public()  (chat-forward public tour)
  page_help.py      ── Help / how-it-works / provenance / limitations
"""
from __future__ import annotations
from pathlib import Path
import sys as _sys
import streamlit as st

APP_DIR = Path(__file__).parent
if str(APP_DIR) not in _sys.path:
    _sys.path.insert(0, str(APP_DIR))

import app_core  # noqa: E402
import page_home  # noqa: E402
import page_researcher  # noqa: E402
import page_discover  # noqa: E402
import page_public  # noqa: E402
import page_help  # noqa: E402


st.set_page_config(
    page_title="AgeFold — Conversational AI for the Biology of Aging",
    page_icon=str(APP_DIR / "assets" / "favicon.png"),
    layout="wide",
    initial_sidebar_state="auto",
)

app_core.inject_theme()

# ---- Pages ----
_home = st.Page(page_home.render, title="Home", icon="🏠", url_path="home", default=True)
_researcher = st.Page(page_researcher.render, title="Researcher", icon="🔬", url_path="researcher")
_discover = st.Page(page_discover.render, title="Discover", icon="🔍", url_path="discover")
_public = st.Page(page_public.render, title="Public", icon="🌍", url_path="public")
_help = st.Page(page_help.render, title="Help", icon="❓", url_path="help")

# Stash Page objects so any page can jump via st.switch_page (Home/Help CTAs).
st.session_state["_nav_pages"] = {
    "home": _home,
    "researcher": _researcher,
    "discover": _discover,
    "public": _public,
    "help": _help,
}

# Native nav is hidden — each page renders its own persistent top-nav row
# (app_core.top_nav) so Home/Researcher/Discover/Public/Help links are always visible
# and you can jump back to Home from any page.
_nav = st.navigation([_home, _researcher, _public, _discover, _help], position="hidden")
_nav.run()
