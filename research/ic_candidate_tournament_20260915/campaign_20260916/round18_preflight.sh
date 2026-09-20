#!/bin/bash
# Wait for round-0018 prepare to finish, check the frozen contract against the
# pre-registration, and only then launch the run. Refuses to run if anything
# the pre-registration fixed came out different.
set -u
cd /home/user/crypto/research/ic_candidate_tournament_20260915 || exit 1

while [ ! -f runs/round-0018b/contract.json ] || [ ! -f runs/round-0018b/fixtures.json ]; do
  if ! pgrep -f "tournament.py prepare --out runs/round-0018b" > /dev/null; then
    echo "PREPARE EXITED WITHOUT A CONTRACT"
    tail -30 runs/round-0018b/build.log
    exit 1
  fi
  sleep 30
done
echo "=== prepare complete at $(date -u +%H:%M:%S)"

python3 - <<'PY'
import json, sys
c = json.load(open('runs/round-0018b/contract.json'))
arms = [a['id'] for a in json.load(open('runs/round-0018b/candidates.json'))]
fx = {k: len(v) for k, v in json.load(open('runs/round-0018b/fixtures.json')).items()}
want = {'cells': ['n13a0','n17a1','n19a0','n23a0','n23a1','n31a0'],
        'holdout_cells': ['n19a1','n29a1'], 'seed': 2026091818,
        'profile': 'pilot', 'target_count': 1, 'objective': 'rho'}
bad = []
for k, v in want.items():
    got = c.get(k)
    print(f'  {k:16} = {got}')
    if got != v:
        bad.append(f'{k}: got {got!r}, pre-registered {v!r}')
print(f'  {"arms":16} = {arms}')
print(f'  {"fixtures":16} = {fx}')
if arms != ['incumbent', 'block', 'column', 'both']:
    bad.append(f'arms: got {arms}')
# The pre-registration budgets 2,340 trials: 36 aa, 90 smoke, 270 development,
# 216 selection, 864 confirmation, 864 replay, with four IC arms plus rho.
planned = 36 + 90 + 270 + 216 + 864 + 864
lim = c.get('limits', {})
cap = lim.get('max_profiled_jobs')
print(f'  pre-registered trial budget = {planned}, max_profiled_jobs = {cap}')
print(f'  {"limits":16} = cpu {lim.get("cpu")}, timeout {lim.get("timeout_seconds")}s')
if cap is None or planned > int(cap):
    bad.append(f'budget {planned} against max_profiled_jobs {cap}')
if lim.get('cpu') != 3 or float(lim.get('timeout_seconds', 0)) != 60.0:
    bad.append(f'limits differ from the pre-registration: {lim}')
if bad:
    print('\nCONTRACT DOES NOT MATCH THE PRE-REGISTRATION:')
    for b in bad:
        print('  -', b)
    sys.exit(2)
print('\ncontract matches the pre-registration')
PY
rc=$?
[ $rc -ne 0 ] && { echo "REFUSING TO RUN (exit $rc)"; exit $rc; }

# THE BASELINE. Round 0018's first attempt was given round 0017's *incumbent*
# instead of round 0017's promoted winner, and every arm-versus-incumbent
# number it produced conflated this round's levers with round 0017's
# certificate. The guard checked eight other fields and not this one, so it
# checks it first now: the incumbent worker must be the round-0017 winner.
WINNER=/home/user/crypto/research/ic_candidate_tournament_20260915/runs/round-0017/source_candidates/orbits/worker
if [ "$(sha256sum < runs/round-0018b/worker | cut -d' ' -f1)" != "$(sha256sum < "$WINNER" | cut -d' ' -f1)" ]; then
  echo "REFUSING TO RUN: the incumbent is not the round-0017 winner (orbits)"
  echo "  built:  $(sha256sum < runs/round-0018b/worker | cut -c1-16)"
  echo "  winner: $(sha256sum < "$WINNER" | cut -c1-16)"
  exit 3
fi
echo "  incumbent: is the round-0017 winner (orbits)"
if ! grep -q factor_base_orbits runs/round-0018b/source/examples/ic_tournament_worker.rs; then
  echo "REFUSING TO RUN: the baseline source does not carry the orbit certificate"
  exit 3
fi
echo "  baseline source carries the orbit-representative certificate"

# The workers the tournament built must be the binaries development measured.
# ROUND18_ARM_BINS, when set, holds the arm workers the development
# measurement ran, so the tournament's own builds can be checked against them.
# Unset is not a failure -- the check is skipped and says so -- because the
# binaries are build outputs and are not committed.
BIN="${ROUND18_ARM_BINS:-}"
if [ -z "$BIN" ]; then
  echo "  (ROUND18_ARM_BINS unset: not comparing built workers to the measured binaries)"
fi
for a in block column both; do
  w="runs/round-0018b/source_candidates/$a/worker"
  m="$BIN/arm18-$a"
  if [ -n "$BIN" ] && [ -f "$w" ] && [ -f "$m" ]; then
    if [ "$(sha256sum < "$w" | cut -d' ' -f1)" = "$(sha256sum < "$m" | cut -d' ' -f1)" ]; then
      echo "  $a worker: identical to the binary development measured"
    else
      echo "  $a worker: DIFFERS from the binary development measured"
    fi
  fi
done

echo "=== launching run at $(date -u +%H:%M:%S)"
exec python3 tournament.py run --round runs/round-0018b --stage all
