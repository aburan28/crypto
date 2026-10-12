#!/bin/bash
set -euo pipefail
# Final builds use the root-owned, image-pinned Cargo configuration.
export CARGO_HOME=/opt/cargo
python3 /opt/harness/final_verify.py
