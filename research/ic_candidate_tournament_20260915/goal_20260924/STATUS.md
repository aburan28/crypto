# September 24 bounded IC goal

Status: implementation/admission work; no new performance tournament launched.
Foundation: [PR 704](https://github.com/aburan28/crypto/pull/704), merged at
`c04879dbb305423d4e7bad9b9e9dd2188f000a35`.

## Objective and stopping rule

Complete at most three improvement rounds. Target at least 20% lower complete
cold instruction cost **and** native cold time than the strongest qualified
archived IC incumbent, on the declared workload only. Finish with a qualifying
winner or an audited finding that none qualified in those rounds. Scaffolding,
an incomplete campaign or exhausted implementation time is not completion.

Before measurement, freeze at least five exact curve/subgroup cells, at least
60 fresh single-target confirmation cases, three process repetitions per target,
resources, accounting unit, floor/reference and stopping rules. Up to 16 complete
pipelines per round; retain six diverse candidates including an exploration
slot. Measure actual combinations, including mechanisms with individually slower
parents. One challenger per round enters confirmation; never retune on it.

Predeclare a familywise 95% rule across at most three confirmation attempts.
Promotion requires both 20% point improvements, uncertainty excluding regression,
no per-cell regression above 10%, independently certified answers for every
target, and passing confirmation and replay. The precise allocation, estimator,
familywise rule and panel must be sealed before the first new measured round;
this status note does not substitute for that protocol.

## Baseline inventory, not yet fresh qualification

| Source | Durable identity | Why retained |
|---|---|---|
| Round 0020 `both` | Source manifest `69de47e30a267e93ab6da91903217062f2e2d3fbe9d6963e17700f77f33c3bd6` | Last formally promoted single-target winner |
| Round 0023 `scaled` | Source manifest `55154f73c35b1f55240b35fd1a5e8e2488c6df41444b5114949e41741c0c3db1` | Improved the archived incumbent; nonpromotion under the old rho objective does not disqualify it as a stronger IC reference |
| Round 0024 `pairinv` proposal | `campaign_20260916/round24-pairinv.patch` over `scaled` | Retained source change and public equivalence check; inspect and qualify before choosing a baseline |

The archive manifest and restore tool are in `../evidence/`. Archive bytes:

- Round 0020 SHA-256 `944a4bd2d6e54fc95fc7bb108ef09a5a1319da45882cf5565c30a8d30f653c2f`.
- Round 0023 SHA-256 `e2aa7111bb03ae606a3af6f249c9bd18beb3df1cc3ca4f813812010336dd47da`.

Both archives restored successfully during this audit. Restoration checks hashes;
it is not a fresh correctness/performance qualification. Preserve the original
sealed sources. Instrument derived copies and measure observer overhead and
equivalence against originals on public development fixtures. Also inspect the
current `icx` engine before claiming the strongest compatible incumbent.

## Next gates

1. Canonical record and base census regression/CI: implemented in `identity.py`,
   `measurement.py`, `test_records.py` and `ci_smoke.py`; review and merge.
2. Instrument the actual optimized producer into eleven exclusive phases,
   retaining ordinary-query outcomes/rank and matrix diagnostics. Combined old
   labels remain unknown under the new schema. Integrate admission into both
   development and promotion drivers before any new comparison.
3. Restore/qualify candidate IC sources and a strong matched rho reference;
   freeze rho width and source before held-out data.
4. Seal the panel and familywise protocol, then run bounded rounds. Keep all
   failures and source/fixture/profiler artifacts; archive them durably, update
   the existing scoreboard, and merge implementation/evidence PRs.

The local development machine is macOS arm64. Calibrated measurements require
the existing Linux amd64 / Valgrind 3.22.0 workflow; native local timings cannot
substitute for that accounting model. No new performance result is claimed here.
