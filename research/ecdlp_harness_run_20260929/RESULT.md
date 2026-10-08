# Result: ECDLP harness run (2026-09-29)

**Verdict:** harness healthy on ECDLP-related problems offline. All five
selected eval tasks passed; `EXP-SEMAEV-001` wrote 12/12 `completed_valid`
records; walkviz recovered `k=118` on the 16-bit seed-7 instance.

**Class:** accounting / operational — no algorithm change in this repository.

## Table (one unit: verified completion)

| variant | verified | notes |
| --- | ---: | --- |
| EVAL-CAP-DLOG-12 (rho) | yes | `k=9`, discrete_log_certificate |
| EVAL-CAP-DLOG-16 (rho) | yes | `k=118`, discrete_log_certificate |
| EVAL-CAP-CURVE-ORDER | yes | order 19 on F₁₇ textbook curve |
| EVAL-CAP-VERIFY-SCRIPT | yes | `k=11`, verify.py exit 0 |
| EVAL-DISC-NO-SOLUTION | yes | honest `solved: false` |
| EXP-SEMAEV-001 rho ×6 | yes | bits 8/10/12 × seeds 1/2 |
| EXP-SEMAEV-001 gb ×6 | yes | 3 decompositions found, 3 trivial ideals |
| walkviz DP + kangaroo | yes | `secret_k_matches: true`, k=118 |

## Doctor / unit tests

- `orchestration doctor`: Ready (inference backends unset; eval agent runtime skipped).
- Suite validate (capability + discipline) on `local`: OK.
- `tests/test_harness.py`: 57 passed, 1 skipped, 1 failed
  (`test_md5_pin_mechanism_real_registry_is_distinct` — MD5 dual-impl probe,
  unrelated to ECDLP).
- This repo's `autolab.py doctor`: reports host Linux/x86_64;
  `native_build_ready=false` without a prepared producer source tree
  (expected for a bare checkout).

## What was not run

- Token-spending `orchestration.eval run` / `make loop` (no API keys).
- Harbor AutoLab improvement campaign (`research/ic_autolab_harness_20260915`)
  — requires the frozen Docker image not present in this environment.
- IC candidate tournament promotion measurement.

## Evidence

Compact receipts under [`evidence/`](evidence/). Full Semaev run directories and
SVG traces were retained in the agent artifact store for this session.

Companion commit: `3e12d32e845a682d685db0a90c0a3847ea159058`.
