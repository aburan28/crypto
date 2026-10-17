#!/bin/bash
# Operational wrapper only: all scheduling/measurement logic is native Rust.
set -u
cd /workspace/repos/crypto-n83-basis || exit 98
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_basis_setup/training.reservation.jsonl --label n83-basis-training-20261015 -- /workspace/n83_basis_setup/frozen_binaries/runner run /workspace/n83_basis_setup/training-plan/matrix.json /workspace/n83_basis_setup/training-campaign > /workspace/n83_basis_setup/training.stdout.jsonl 2> /workspace/n83_basis_setup/training.stderr.log
n83_training_exit=$?
( set -o noclobber; printf '%s\n' "$n83_training_exit" > /workspace/n83_basis_setup/training.exit.txt )
exit "$n83_training_exit"
