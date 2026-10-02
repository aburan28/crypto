# Result: ECDLP harness runs (2026-10-01)

**Verdict:** harness healthy across seven consecutive offline executions.
All 7 × 12 EXP-SEMAEV-001 records are `completed_valid`; walkviz closes
identically each time; `tests/test_harness.py` passes 59/59 every run.

**Class:** accounting / operational — no algorithm change in this repository.

## Table (one unit: verified completion)

| execution | companion commit (executed) | tree | EXP-SEMAEV-001 | walkviz | pytest |
| --- | --- | --- | --- | --- | --- |
| run 1 | `4129b1791d25` | dirty | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |
| run 2 | `4129b1791d25` | dirty | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |
| run 3 | `4129b1791d25` | dirty | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |
| run 4 | `4129b1791d25` | dirty | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |
| run 5 | `4129b1791d25` | dirty | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |
| run 6 | `7b6017bd7df3` | dirty | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |
| run 7 | `fac5b119772c` | clean | 12/12 valid | closed, 31 ops / tail 28 | 59 passed |

Per-record composition (inspected in full on run 1; runs 2–7 checked at
manifest-status level, all `completed_valid`): 6 rho discrete-log
certificates verified by independent recompute; 3 GB decompositions
verified (b8-s1, b8-s2, b10-s2); 3 trivial-ideal no-claim records verified.

The companion checkout advanced externally mid-session (4129b17 → 7b6017b →
fac5b11); results are identical across all three commits. Doctor + suite
validation (`capability` 4/4, `discipline` 8/8, local backend) were executed
on the first two sessions; experiment/walkviz/pytest on all seven.

## Boundaries

- Toy-scale verification only. No ECC2K-130, m=83, frozen-tournament, or
  end-to-end IC `S` claim.
- No token-spending `orchestration.eval run` / agent loop (offline path).
- Evidence dirs are local-only; content hashes and byte counts are recorded
  in [`evidence/manifest.json`](evidence/manifest.json).

Companion repo: `aburan28/crypto-autoresearcher`
(`4129b1791d252f40942337c8f66ecb7484b55e91`,
`7b6017bd7df3`, `fac5b11977` — per-execution commit above).
