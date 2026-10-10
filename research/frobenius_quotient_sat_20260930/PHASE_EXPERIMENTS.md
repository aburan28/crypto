# Phase-aware representation experiment plan

## Existing verified baseline
The existing `README.md`, `run.py`, `encode.py`, `test_encode.py`, `brute_check.py`, and `summarize.py` remain the source of truth. The published-in-repo trials show n=13 verified solutions with unknown phases but n=19 unknown-phase planted cases timing out at 1200 s. These are stage diagnostics, not ECDLP speedups.

## Added experiment runner
```sh
cd research/frobenius_quotient_sat_20260930
pip install pycryptosat==5.16.0
python3 -m unittest test_encode -v
python3 phase_experiments.py
python3 phase_experiments.py --execute --fields 13 --targets 2 --budget 30 --procs 1
python3 summarize.py phase_results/n13_phase_subgroup_xor.jsonl
python3 brute_check.py phase_results/n13_phase_subgroup_xor.jsonl
```
The first invocation is a dry run. All executions record a receipt; each run reuses the original exact-verification solver. Timeouts are expected and counted. The matrix is *not* a statistically controlled comparison of identical target instances between phase-on and phase-off; phase-on planted target generation uses randomized Frobenius shifts. For causal comparisons, add a frozen target manifest and use identical target coordinates.

## Next encoding experiments (not yet implemented)
1. **Binary phase selectors:** represent k in ceil(log2 n) bits, constrain k < n, and implement controlled Frobenius powers with Boolean gates. Compare actual CNF/XOR size and end-to-end cost, not variable count alone.
2. **Relative phases:** exploit common Frobenius action only when the target equation is simultaneously transformed, or target is fixed by the action. Do not arbitrarily set one phase to zero for a fixed arbitrary target.
3. **Orbit canonicalization:** encode representative and phase metadata with a round-trip identity check; handle short orbits and sign carefully.
4. **Matrix/rank integration:** independently verify each signed point relation and the coefficient weights induced by Frobenius eigenvalues modulo subgroup order.
5. **Scaling:** move to n=19 only after an n=13 improvement reproduces across many seeds. n=29 and n=83 are gated on n=19 success.

## Metrics and acceptance
Record build/load/solve/verification CPU and wall time, timeouts, UNSAT, rejected tuples, gates, variables, clauses, verified relations, duplicate relations, rank, and full precomputation. Compare identical target sets, independent seeds, interleaved order, and confidence intervals. A paper claim requires a reproducible end-to-end gain; no speedup is claimed here.
