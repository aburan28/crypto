#!/bin/bash
# Wait for round-0019 prepare to finish, check the frozen contract against the
# pre-registration, and only then launch the run. Refuses to run if anything
# ROUND19-single-target.md fixed came out different.
set -u
cd /home/user/crypto/research/ic_candidate_tournament_20260915 || exit 1
ROUND=runs/round-0019

while [ ! -f $ROUND/contract.json ] || [ ! -f $ROUND/fixtures.json ]; do
  if ! pgrep -f "tournament.py prepare --out $ROUND" > /dev/null; then
    echo "PREPARE EXITED WITHOUT A CONTRACT"
    tail -30 $ROUND/build.log 2>/dev/null
    exit 1
  fi
  sleep 30
done
echo "=== prepare complete at $(date -u +%H:%M:%S)"

python3 - <<'PY'
import json, sys, collections
ROUND = 'runs/round-0019'
c = json.load(open(f'{ROUND}/contract.json'))
arms = [a['id'] for a in json.load(open(f'{ROUND}/candidates.json'))]
fixtures = json.load(open(f'{ROUND}/fixtures.json'))
fx = {k: len(v) for k, v in fixtures.items()}
want = {'cells': ['n13a0','n17a1','n19a0','n23a0','n23a1','n31a0'],
        'holdout_cells': ['n19a1','n29a1'], 'seed': 2026092119,
        'profile': 'pilot', 'target_count': 1, 'objective': 'rho',
        'confirmation_cases': 124}
bad = []
for k, v in want.items():
    got = c.get(k)
    print(f'  {k:26} = {got}')
    if got != v:
        bad.append(f'{k}: got {got!r}, pre-registered {v!r}')
print(f'  {"arms":26} = {arms}')
print(f'  {"fixtures":26} = {fx}')
if arms != ['incumbent', 'both']:
    bad.append(f'arms: got {arms}, pre-registered [incumbent, both]')

# THE ALLOCATION. This is the round's one protocol change, so it is checked
# twice: against the contract's own record of it, and against the fixtures the
# contract produced, which is what the run will actually measure.
alloc = {'n13a0': 12, 'n17a1': 12, 'n19a0': 12, 'n19a1': 12,
         'n23a0': 12, 'n23a1': 40, 'n29a1': 12, 'n31a0': 12}
got = c.get('confirmation_cases_per_cell')
print(f'  {"allocation (contract)":26} = {got}')
if got != alloc:
    bad.append(f'allocation: got {got!r}, pre-registered {alloc!r}')
for stage in ('confirmation', 'replay'):
    drawn = collections.Counter(case['cell'] for case in fixtures[stage])
    print(f'  {"allocation (" + stage + ")":26} = {dict(sorted(drawn.items()))}')
    if dict(drawn) != alloc:
        bad.append(f'{stage} fixtures: got {dict(sorted(drawn.items()))!r}, pre-registered {alloc!r}')
    if len({case['fixture_sha256'] for case in fixtures[stage]}) != sum(alloc.values()):
        bad.append(f'{stage}: repeated fixtures, which would understate the spread')

# The pre-registration budgets 2,646 trials: 36 aa, 54 smoke, 162 development,
# 162 selection, 1,116 confirmation, 1,116 replay, with two IC arms plus rho.
planned = 36 + 54 + 162 + 162 + 1116 + 1116
lim = c.get('limits', {})
cap = lim.get('max_profiled_jobs')
print(f'  pre-registered trial budget = {planned}, max_profiled_jobs = {cap}')
print(f'  {"limits":26} = cpu {lim.get("cpu")}, timeout {lim.get("timeout_seconds")}s')
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
# number it produced conflated that round's levers with round 0017's
# certificate. The guard checked eight other fields and not this one.
WINNER=runs/round-0017/source_candidates/orbits/worker
if [ -f "$WINNER" ]; then
  if [ "$(sha256sum < $ROUND/worker | cut -d' ' -f1)" != "$(sha256sum < "$WINNER" | cut -d' ' -f1)" ]; then
    echo "REFUSING TO RUN: the incumbent is not the round-0017 winner (orbits)"
    echo "  built:  $(sha256sum < $ROUND/worker | cut -c1-16)"
    echo "  winner: $(sha256sum < "$WINNER" | cut -c1-16)"
    exit 3
  fi
  echo "  incumbent: is the round-0017 winner (orbits), retained by round 0018b"
else
  echo "  (round-0017 winner not restored: comparing the incumbent to round 0018b instead)"
  if [ "$(sha256sum < $ROUND/worker | cut -d' ' -f1)" != "$(sha256sum < runs/round-0018b/worker | cut -d' ' -f1)" ]; then
    echo "REFUSING TO RUN: the incumbent is not the round-0018b retained winner"
    exit 3
  fi
fi
if ! grep -q factor_base_orbits $ROUND/source/examples/ic_tournament_worker.rs; then
  echo "REFUSING TO RUN: the baseline source does not carry the orbit certificate"
  exit 3
fi
echo "  baseline source carries the orbit-representative certificate"

# The challenger must carry both round-0018 patches and declare the convention
# the amended oracle reads; an arm that silently lost a patch would read as a
# null result rather than as a broken build.
CH=$ROUND/source_candidates/both/source/examples/ic_tournament_worker.rs
if [ -f "$CH" ]; then
  grep -q 'column_convention.*representative' "$CH" || {
    echo "REFUSING TO RUN: the challenger does not declare column_convention: representative"; exit 3; }
  grep -q 'clamp(8, 16)' $ROUND/source_candidates/both/source/src/cryptanalysis/koblitz_tiny_ic.rs || {
    echo "REFUSING TO RUN: the challenger does not carry the scan-block ceiling"; exit 3; }
  echo "  challenger carries both round-0018 patches"
fi

echo "=== launching run at $(date -u +%H:%M:%S)"
exec python3 tournament.py run --round $ROUND --stage all
