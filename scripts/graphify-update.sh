#!/usr/bin/env bash
# Rebuild scoped graphify graph (mofography/ only) -> graphify-out/
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
rm -rf mofography/graphify-out
graphify update mofography --no-viz
rm -rf graphify-out
mv mofography/graphify-out graphify-out
NODES=$(python3 -c "import json; print(len(json.load(open('graphify-out/graph.json'))['nodes']))")
echo "graphify-out: ${NODES} nodes (mofography/ scope)"
