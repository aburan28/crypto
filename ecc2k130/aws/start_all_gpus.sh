#!/usr/bin/env bash
#
# One campaign worker per visible GPU. EC2 bootstrap already starts
# ecc2k130-worker@N through systemd; this is the RunPod / no-systemd path.
#
# An 8-wide B200 MIG pod that runs a single process shows ~1/8 GPU
# utilization. This fills every device nvidia-smi lists.
#
#   ECC_ROOT=/opt/ecc2k130 ./aws/start_all_gpus.sh
#
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
export ECC_ROOT=${ECC_ROOT:-$(cd "$HERE/.." && pwd)}
export PYTHONUNBUFFERED=1
exec env ECC_ALL_GPUS=1 python3 -u "$HERE/worker.py"
