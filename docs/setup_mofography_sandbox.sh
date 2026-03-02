#!/usr/bin/env bash

set -euo pipefail

# Simple sandbox setup script for the mofography project using uv.
# Usage:
#   ./setup_mofography_sandbox.sh [TARGET_DIR]
# If TARGET_DIR is not provided, it defaults to "mofography".
#
# Authentication:
#   For private GitHub access, set the GITHUB_PAT environment variable
#   to a token that can read the arung-northwestern/mofography repo.
#   If GITHUB_PAT is not set, the script will prompt for it.

TARGET_DIR="${1:-mofography}"

ROOT_DIR="$(pwd)"
ABS_TARGET_DIR="${ROOT_DIR}/${TARGET_DIR}"

echo "Creating mofography sandbox in '${ABS_TARGET_DIR}'"

if ! command -v python3 >/dev/null 2>&1; then
  echo "Error: python3 is required but was not found on PATH." >&2
  exit 1
fi

if ! command -v git >/dev/null 2>&1; then
  echo "Error: git is required but was not found on PATH." >&2
  exit 1
fi

if ! command -v cargo >/dev/null 2>&1; then
  echo "Error: Rust toolchain (cargo) is required to build the mgr_rust backend." >&2
  echo "Please install Rust from https://rustup.rs/ and re-run this script." >&2
  exit 1
fi

if ! command -v uv >/dev/null 2>&1; then
  echo "'uv' was not found on PATH; attempting to install it with pip..."
  python3 -m pip install --user --upgrade uv
  # Ensure the typical user base install directory is on PATH
  export PATH="$HOME/.local/bin:$PATH"
fi

if ! command -v uv >/dev/null 2>&1; then
  echo "Error: 'uv' is still not available after installation attempt." >&2
  echo "Please install uv manually (see https://github.com/astral-sh/uv) and re-run this script." >&2
  exit 1
fi

mkdir -p "${ABS_TARGET_DIR}"
cd "${ABS_TARGET_DIR}"

echo "Creating virtual environment with uv in '.venv'..."
uv venv .venv

# shellcheck disable=SC1091
source .venv/bin/activate

echo "Upgrading pip and installing maturin..."
uv pip install --upgrade pip
uv pip install maturin

REPO_OWNER="arung-northwestern"
REPO_NAME="mofography"
BASE_REPO_URL="https://github.com/${REPO_OWNER}/${REPO_NAME}.git"

GITHUB_PAT_ENV="${GITHUB_PAT:-}"

if [ -z "${GITHUB_PAT_ENV}" ]; then
  echo "A GitHub personal access token with read access to ${REPO_OWNER}/${REPO_NAME} is required."
  read -r -s -p "Enter GitHub personal access token: " GITHUB_PAT_INPUT
  echo
  GITHUB_PAT_ENV="${GITHUB_PAT_INPUT}"
fi

if [ -z "${GITHUB_PAT_ENV}" ]; then
  echo "Error: GitHub personal access token was not provided." >&2
  exit 1
fi

CLONE_URL="https://${GITHUB_PAT_ENV}@github.com/${REPO_OWNER}/${REPO_NAME}.git"

if [ -d "mofography-src/.git" ]; then
  echo "Directory 'mofography-src' already exists and appears to be a git repository; skipping clone."
else
  echo "Cloning ${REPO_OWNER}/${REPO_NAME} into 'mofography-src'..."
  git clone "${CLONE_URL}" mofography-src
fi

# Clear the token from variables as best we can
unset GITHUB_PAT_ENV
unset GITHUB_PAT_INPUT

echo "Building mgr_rust backend with maturin (release)..."
cd "${ABS_TARGET_DIR}/mofography-src/mgr_rust"
maturin build --release

WHEEL_PATH="$(ls -1 target/wheels/mgr_rust-*.whl 2>/dev/null | head -n 1 || true)"
if [ -z "${WHEEL_PATH}" ] || [ ! -f "${WHEEL_PATH}" ]; then
  echo "Error: mgr_rust wheel was not found under 'target/wheels/'." >&2
  exit 1
fi

echo "Installing mgr_rust wheel into the virtual environment..."
uv pip install "${WHEEL_PATH}"

echo "Installing Python package 'mofography' from source..."
cd "${ABS_TARGET_DIR}/mofography-src"
uv pip install .

echo "Installing additional Python dependencies with uv..."
uv pip install \
  pymatgen \
  ase \
  dask \
  dask-jobqueue \
  networkx \
  scikit-image \
  tqdm \
  aiohttp \
  sparse \
  scikit-learn \
  rich \
  trimesh \
  pyvista \
  plotly \
  furl \
  requests \
  freud-analysis \
  kaleido \
  rdkit \
  optuna \
  stackprinter

echo
echo "==============================================================="
echo "Mofography sandbox setup complete."
echo
echo "To activate the virtual environment in a new shell, run:"
echo "  source \"${ABS_TARGET_DIR}/.venv/bin/activate\""
echo
echo "Then you can run the reproduction script, for example:"
echo "  cd \"${ROOT_DIR}\""
echo "  source \"${ABS_TARGET_DIR}/.venv/bin/activate\""
echo "  python repro_segmentation_results.py"
echo "==============================================================="

