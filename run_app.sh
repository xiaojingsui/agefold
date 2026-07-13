#!/usr/bin/env bash
# Launch the LiP-MS Aging Viewer locally on your Mac.
# Uses the conda env at ~/.claude-science/conda/envs/lipms-viewer directly
# (no `conda activate` needed).
#
# Prerequisites:
#   1. ANTHROPIC_API_KEY in .streamlit/secrets.toml (or in your shell env)
#   2. AlphaFold PDBs on OneDrive (referenced by app/data/alphafold_index.parquet)

set -euo pipefail
cd "$(dirname "$0")"

STREAMLIT_BIN="$HOME/.claude-science/conda/envs/lipms-viewer/bin/streamlit"

if [ ! -x "$STREAMLIT_BIN" ]; then
  echo "❌ Cannot find $STREAMLIT_BIN"
  echo "   The lipms-viewer conda env is missing. Ask Claude Science to rebuild it."
  exit 1
fi

# API-key check (warn only — offline fallback still works)
if [ -z "${ANTHROPIC_API_KEY:-}" ] && [ ! -f .streamlit/secrets.toml ]; then
  echo "⚠️  No API key found — chat will use the offline report-card fallback."
  echo "   Paste your key into .streamlit/secrets.toml to enable the chat."
elif [ -f .streamlit/secrets.toml ] && grep -q "sk-ant-YOUR-KEY-HERE" .streamlit/secrets.toml; then
  echo "⚠️  .streamlit/secrets.toml still has the placeholder key — chat will error."
  echo "   Open it and paste your real key:"
  echo "     open -e .streamlit/secrets.toml"
fi

echo "▶ launching Streamlit on http://localhost:8501"
echo "▶ ctrl-C to stop"
exec "$STREAMLIT_BIN" run app/streamlit_app.py \
  --server.port 8501 \
  --server.address localhost \
  --browser.gatherUsageStats false
