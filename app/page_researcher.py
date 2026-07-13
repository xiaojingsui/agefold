"""
page_researcher.py — the "Researcher" page.

The full technical viewer: sidebar protein/condition/painting controls plus the
four tabs (Analysis multi-agent workspace, Overview 3D+heatmap+tracks+card+
peptides+ortholog panel, Ask-about-this-protein chat, Discover). Thin wrapper
that hands the shared ctx to researcher_view.render().
"""
from __future__ import annotations
import app_core
import researcher_view


def render() -> None:
    app_core.top_nav("researcher")
    ctx = app_core.build_ctx()
    researcher_view.render(ctx)
