#!/usr/bin/env bash
# Build or incrementally update code-review-graph for blast-radius queries.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p .crg-out
if code-review-graph status --repo "$ROOT" --json 2>/dev/null | python3 -c "import sys,json; d=json.load(sys.stdin); exit(0 if d.get('nodes') else 1)"; then
  code-review-graph update --repo "$ROOT" -q
else
  code-review-graph build --repo "$ROOT" -q
fi
code-review-graph status --repo "$ROOT" --json
