# Run report: run-20261011T1955Z-verify

A descriptive record of the second verification batch for the Eigen-Séta v2
Rust KEM. It states observations only. Interpretation is in `../../README.md`,
`eigen_seta/SIDE_CHANNELS.md` and the paper in `papers/eigen-seta/`. The PDF
is `report.pdf`.

## What ran

| Field | Value |
|---|---|
| Run id | `run-20261011T1955Z-verify` |
| Start, end of phase 1, end (UTC) | 2026-10-11T19:56:36Z, 20:29:49Z, 20:37:48Z (`records/start.txt`, `records/phase1_end.txt`, `records/end.txt`) |
| Sources | Commit `d504af28e` (`records/commit.txt`), clean worktree for `eigen_seta/` (`records/worktree_status.txt` is empty) |
| Build | rustc 1.97.0, release profile with debuginfo (`records/build.log`) |
| Tools | valgrind 3.22 (callgrind), glibc malloc |
| Host | Intel Xeon Processor @ 2.80 GHz, 4 logical cores, Linux 6.18.44 x86_64, shared cloud VM |
| Script | `records/verify_run.sh` |

Phase 1 ran two things in parallel:
* four callgrind probes at L1;
* a native failure-rate stream: `run Toy16 3000`, then `run Toy32 600`, then
  `run L1 20`.

Phase 2 ran the 12-run Toy16 probe grid in three streams.

Phase 3 ran the timing test `dudect_decaps 1000`, alone on the machine.

The probe command was:

```
valgrind --tool=callgrind --toggle-collect='*probe_core*' ./ct_probe <key> <ct> <valid|foreign> <set>
```

Only `decaps_core` runs inside `probe_core`. A foreign ciphertext is well
formed but made for another key, so the FO check rejects it.
`records/ctsummary.sh` splits each count into allocator and the rest.
Allocator means glibc malloc and arena code, heap growth (`sbrk`, `brk`,
`morecore`) and the `mem*` routines.

**Filter-rate step, run right after this batch.** It used
`examples/f_filter_rate.rs` at commit `72051ac67`, which differs from
`d504af28e` only by capping the public basis search and adding that example.
It ran from 2026-10-11T20:38:30Z to 20:46:52Z, with the timestamps inside
`records/f_filter_rate.txt`.

Requested and completed:

| Step | Requested | Completed | Failures |
|---|---|---|---|
| L1 probes | 4 | 4 | none |
| Toy16 probes | 12 | 12 | none |
| Failure-rate trials | 3,620 | 3,620 | none |
| Timing test | 1 | 1 | none |
| Filter-rate measurements | 5 | 5 | none |

## Measured values

### Decapsulation correctness

Records are in `records/failrate/`.

| Set | Trials | Correct | Failures |
|---|---|---|---|
| Toy16 | 3,000 | 3,000 | 0 |
| Toy32 | 600 | 600 | 0 |
| L1 | 20 | 20 | 0 |

### Instructions in `decaps_core`

Records are in `records/{toyprobe,l1probe}/ctprobe_summary.txt`, with raw
outputs in `callgrind/*.xz`.

| Set | Key seed | Ciphertext | Accepted | Total | Allocator | Rest |
|---|---|---|---|---|---|---|
| Toy16 | 1 | valid 1 | yes | 1,147,950,979 | 31,741,763 | 1,116,209,216 |
| Toy16 | 1 | valid 2 | yes | 1,147,950,979 | 31,741,763 | 1,116,209,216 |
| Toy16 | 1 | foreign | no | 1,148,059,231 | 31,850,015 | 1,116,209,216 |
| Toy16 | 2 | valid 1 | yes | 1,147,388,400 | 31,179,184 | 1,116,209,216 |
| Toy16 | 2 | valid 2 | yes | 1,147,388,400 | 31,179,184 | 1,116,209,216 |
| Toy16 | 2 | foreign | no | 1,147,693,051 | 31,483,835 | 1,116,209,216 |
| Toy16 | 3 | valid 1 | yes | 1,147,459,087 | 31,249,871 | 1,116,209,216 |
| Toy16 | 3 | valid 2 | yes | 1,147,459,087 | 31,249,871 | 1,116,209,216 |
| Toy16 | 3 | foreign | no | 1,147,849,768 | 31,640,552 | 1,116,209,216 |
| Toy16 | 4 | valid 1 | yes | 1,147,477,368 | 31,268,152 | 1,116,209,216 |
| Toy16 | 4 | valid 2 | yes | 1,147,477,368 | 31,268,152 | 1,116,209,216 |
| Toy16 | 4 | foreign | no | 1,147,787,126 | 31,577,910 | 1,116,209,216 |
| L1 | 1 | valid 1 | yes | 124,250,795,358 | 3,254,365,748 | 120,996,429,610 |
| L1 | 1 | valid 2 | yes | 124,250,795,358 | 3,254,365,748 | 120,996,429,610 |
| L1 | 1 | foreign | no | 124,257,288,150 | 3,260,858,540 | 120,996,429,610 |
| L1 | 2 | valid 1 | yes | 124,262,712,389 | 3,266,282,779 | 120,996,429,610 |

Observations:
* At each parameter set, the rest takes a single value.
* The allocator count varies with the key and the mode, and it is equal for
  the two valid ciphertexts of each key.
* A first L1 summary did not count `sbrk` and `brk` as allocator code. It
  showed a 174-instruction difference between the keys, all of it in those
  two libc heap-growth routines. The summary script was corrected before
  the final summaries were written.

### Timing test

The record is `records/dudect_output.txt`, in ns per call with 1,000 samples
per class.

```
samples/class 1000; median ns: fixed 122373850, random 120844155
t(raw) = +1.17
t(crop <= p50) = +0.14  (n = 473, 528)
t(crop <= p75) = +1.46  (n = 730, 771)
t(crop <= p90) = +2.36  (n = 899, 902)
```

### Filter acceptance rates

The record is `records/f_filter_rate.txt`. It gives the fraction of odd f
candidates that pass key generation's public filters.

| Set | Drawn | Accepted | Rate | log2 |
|---|---|---|---|---|
| Toy16 | 12,408 | 300 | 0.02418 | −5.37 |
| Toy32 | 24,956 | 300 | 0.01202 | −6.38 |
| L1 | 55,933 | 200 | 0.00358 | −8.13 |
| L3 | 39,987 | 100 | 0.00250 | −8.64 |
| L5 | 51,887 | 100 | 0.00193 | −9.02 |

Most rejections come from gcd(f, 21) ≠ 1 (about 43%) or from the
two-squares step (52–56% of all drawn candidates). The schedule condition rejects roughly half of
the candidates that reach it: 188 of 388 at L1.

## Graphs

Both graphs are drawn from the records by `../../figures/plot_verify.py`.

* `ct_probe.svg` shows, per run, the non-allocator count and the allocator
  count, for Toy16 and L1.
* `f_filter.svg` shows log2 of the filter acceptance rate per set. The error
  bars are the approximate binomial standard error 1/(ln 2 · √k), where k is
  the number of accepted candidates.

## Scope

* The probes cover Toy16 and L1, not L3 or L5.
* One compiler and one target.
* Instruction counts only, not memory addresses.
* The timing test uses one key.
* Zero failures in 3,000 trials bounds the Toy16 failure rate only by about
  10^-3 at 95% confidence.
