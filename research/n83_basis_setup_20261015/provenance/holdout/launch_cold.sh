#!/bin/bash
set -u
cd /workspace/repos/crypto-n83-basis || exit 98
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_basis_setup/cold.reservation.jsonl --label n83-basis-fresh-cold -- /workspace/n83_basis_setup/frozen_binaries/holdout_runner run /workspace/n83_basis_setup/holdout-plan/matrix.json /workspace/n83_basis_setup/cold-campaign cold > /workspace/n83_basis_setup/cold.stdout.jsonl 2> /workspace/n83_basis_setup/cold.stderr.log
n83_cold_exit=$?
( set -o noclobber; printf '%s\n' "$n83_cold_exit" > /workspace/n83_basis_setup/cold.exit.txt )
exit "$n83_cold_exit"
