"""
page_public.py — the "Public / Explore aging" page.

Chat-forward plain-language tour for a lay or education audience. Thin wrapper
that hands the shared ctx to public_view.render_public(), which routes between
the chat-forward landing (protein search + featured stories + glossary) and a
per-protein plain-language page with its own conversational chat.
"""
from __future__ import annotations
import app_core
import public_view


def render() -> None:
    app_core.top_nav("public")
    ctx = app_core.build_ctx()
    public_view.render_public(ctx)
