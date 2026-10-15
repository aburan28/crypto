#!/bin/bash
set -u
cd /workspace/repos/crypto-n83-basis || exit 98
/workspace/n83_evaluation/frozen_binaries/isolated_bench reserve --cpus 3 --out /workspace/n83_basis_setup/cohort.reservation.jsonl --label n83-basis-fresh-cohort -- /workspace/n83_basis_setup/frozen_binaries/holdout_runner run /workspace/n83_basis_setup/holdout-plan/matrix.json /workspace/n83_basis_setup/cohort-campaign cohort > /workspace/n83_basis_setup/cohort.stdout.jsonl 2> /workspace/n83_basis_setup/cohort.stderr.log
n83_cohort_exit=$?
( set -o noclobber; printf '%s\n' "$n83_cohort_exit" > /workspace/n83_basis_setup/cohort.exit.txt )
exit "$n83_cohort_exit"
