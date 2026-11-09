#!/usr/bin/env bash
# Install mofography + mgr_rust (release build) on Linux/macOS.
#
# One-liners (after this file is on GitHub):
#   # Fresh uv sandbox next to your shell cwd:
#   curl -fsSL https://raw.githubusercontent.com/mofography/mofography/master/scripts/install_mofography.sh | bash -s -- --mode uv-sandbox
#
#   # Use active conda/venv (install deps with pip into current Python):
#   curl -fsSL https://raw.githubusercontent.com/mofography/mofography/master/scripts/install_mofography.sh | bash -s -- --mode pip-existing
#
# From a clone:
#   ./scripts/install_mofography.sh --mode uv-sandbox
#   ./scripts/install_mofography.sh --mode pip-existing --repo-dir .

set -euo pipefail

MODE="uv-sandbox"
REPO_DIR=""
REPO_URL="https://github.com/mofography/mofography.git"
REF="master"
PYTHON_BIN="${PYTHON_BIN:-python3}"

usage() {
  sed -n '2,12p' "$0" | sed 's/^# //'
  echo "Options:"
  echo "  --mode uv-sandbox|pip-existing   (default: uv-sandbox)"
  echo "  --repo-dir PATH                  existing clone; skip git clone"
  echo "  --ref TAG_OR_BRANCH              git checkout after clone (e.g. 2026.10.09.1)"
  echo "  --repo-url URL                   (default: mofography/mofography)"
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --repo-dir) REPO_DIR="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --repo-url) REPO_URL="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    *) echo "Unknown option: $1" >&2; usage 1 ;;
  esac
done

MOFO_DEPS=(
  pymatgen ase dask dask-jobqueue networkx scikit-image tqdm aiohttp
  sparse scikit-learn rich trimesh pyvista plotly furl requests freud-analysis
  rdkit optuna kaleido nbformat icecream
)

need_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Error: '$1' not found on PATH." >&2
    exit 1
  fi
}

need_cmd git
need_cmd cargo

if [[ -z "${REPO_DIR}" ]]; then
  REPO_DIR="${MOFOGRAPHY_INSTALL_DIR:-./mofography}"
  if [[ -d "${REPO_DIR}/.git" ]]; then
    echo "Using existing clone at ${REPO_DIR}"
    cd "${REPO_DIR}"
    git fetch --tags origin 2>/dev/null || true
    git checkout "${REF}"
  else
    echo "Cloning ${REPO_URL} into ${REPO_DIR} (ref: ${REF})..."
    git clone --depth 1 --branch "${REF}" "${REPO_URL}" "${REPO_DIR}" 2>/dev/null \
      || { git clone "${REPO_URL}" "${REPO_DIR}" && cd "${REPO_DIR}" && git checkout "${REF}"; }
    cd "${REPO_DIR}"
  fi
else
  REPO_DIR="$(cd "${REPO_DIR}" && pwd)"
  cd "${REPO_DIR}"
fi

install_python_stack() {
  local pip_cmd=("$@")
  echo "Installing Python dependencies..."
  "${pip_cmd[@]}" install --upgrade pip maturin
  "${pip_cmd[@]}" install "${MOFO_DEPS[@]}"
}

case "${MODE}" in
  uv-sandbox)
    if ! command -v uv >/dev/null 2>&1; then
      echo "Installing uv..."
      need_cmd "${PYTHON_BIN}"
      "${PYTHON_BIN}" -m pip install --user uv
      export PATH="${HOME}/.local/bin:${PATH}"
    fi
    need_cmd uv
    uv venv .venv
    # shellcheck disable=SC1091
    source .venv/bin/activate
    install_python_stack uv pip
    ;;
  pip-existing)
    need_cmd "${PYTHON_BIN}"
    install_python_stack "${PYTHON_BIN}" -m pip
    ;;
  *)
    echo "Unknown --mode: ${MODE}" >&2
    usage 1
    ;;
esac

echo "Building mgr_rust (maturin develop --release)..."
cd mgr_rust
maturin develop --release
cd ..

echo "Installing mofography (editable)..."
python -m pip install -e .

echo ""
echo "Done. Verify:"
echo "  python -c \"import mgr_rust, mofography as m; print(m.__version__)\""
if [[ "${MODE}" == "uv-sandbox" ]]; then
  echo "Activate later: source \"${REPO_DIR}/.venv/bin/activate\""
fi
