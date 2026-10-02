#!/usr/bin/env bash
# Self-check for the TPU index-calculus backend.
#
# Runs the whole suite on CPU (JAX CPU + Pallas interpret mode), so it needs
# no TPU and no device.  This is the host-verification gate: the kernels are
# correct iff they agree with the scalar oracle in ic/reference.py.  It makes
# NO performance claim -- see protocol/RESEARCH_TPU_IC.md.
set -euo pipefail
cd "$(dirname "$0")"

PY="${PYTHON:-python3}"
export JAX_PLATFORMS="${JAX_PLATFORMS:-cpu}"

"$PY" -c "import jax, numpy" 2>/dev/null || {
  echo "jax/numpy not importable; install with: $PY -m pip install -r requirements.txt" >&2
  exit 2
}

exec "$PY" -m pytest -q tests "$@"
