# One-shot frozen-Q point-sum cold panel: diagnostic, no admitted timing

**Decision.** The optional full-point index is not promoted. The one-shot
[workflow run 36764654520](https://github.com/aburan28/crypto/actions/runs/36764654520)
at merged harness commit `53a0238f3e829f1a5c60ec70ce036308241de5fe`
completed all 180 cold child processes, but zero of six cells satisfies the
frozen verifier's full timing gate. This is **not** a formal rho crossover
or no-go result. Method-level speedup, common operation-unit cost and
ECC2K-130 extrapolation remain unset. The prior scoreboard verdict is
unchanged.

The table reports whole-process `wait4` child user+system CPU only as a
**diagnostic**. Each point/control ratio uses the geometric mean of the
two S3 controls within a rotated block. Brackets are paired two-sided 95%
log-t intervals; repeated Q across blocks measures run noise, not Q
distribution uncertainty. Lower than 1 favors point-sum. The final column
is `1-1/(median point/rho)` and is a diagnostic estimate of the remaining
point-sum CPU reduction needed to match rho at that cell.

| Cell | K | Blocks | Control B/A | Point/control | Point/rho | Required reduction | Frozen replay |
|:--|--:|--:|:--|:--|:--|--:|:--|
| n37/L1 | 7 | 20 | 1.003 [0.988, 1.016] | 0.994 [0.990, 1.006] | 2.926 [2.865, 2.966] | 65.8% | PASS |
| n37/L1024 | 42 | 5 | 0.997 [0.988, 1.008] | 1.009 [0.999, 1.020] | 3.243 [3.227, 3.258] | 69.2% | **FAIL** |
| n41/L1 | 85 | 5 | 0.998 [0.960, 1.023] | 0.977 [0.963, 1.003] | 39.129 [36.764, 40.320] | 97.4% | PASS |
| n41/L1024 | 255 | 5 | 1.007 [0.996, 1.016] | 1.095 [1.081, 1.107] | 1.470 [1.454, 1.482] | 32.0% | PASS |
| n53/L1 | 220 | 5 | 1.002 [0.984, 1.017] | 0.969 [0.932, 0.994] | 29.181 [28.650, 29.895] | 96.6% | PASS |
| n53/L1024 | 440 | 5 | 1.001 [0.921, 1.083] | 1.017 [0.966, 1.064] | 1.476 [1.412, 1.542] | 32.3% | PASS |

All six A/A medians are within [0.9, 1.1] and their intervals include 1.
Each reserved-core monitor recorded **zero contended samples**. The
pre-dispatch verifier also required an empty `left_on_reserved.user_threads`
list, a stricter condition than the protocol's sampled-contention wording.
The hosted VMs retained 81–83 user threads with affinity to the reserved
CPUs despite the monitor moving 43–48 other threads. The frozen verifier
therefore marks each otherwise passing cell `timing_eligible=false`.
We retain that stricter pre-dispatch gate for this run; changing the gate
after seeing ratios would change the admission rule retrospectively.
The n53/L1 point/control interval below 1 is likewise only a diagnostic.

The n37/L1024 producer itself finished all 20 children and recovered all
1,024 target logs per child. Its hosted and Mac frozen verifiers both fail
the same `pinned_intermediates` assertion at target index 763. The
[separate diagnostic](N37_BATCH_DIAGNOSTIC.json) independently replays all
15 compact rank traces, 15,360 compact target logs and 5,120 rho logs.
For that one Q, the first and third factor-base x-codes coincide. The
solver's selected lift and a second valid lift both sum to the public Q;
the second lift has the pinned pair x-values. The frozen verifier binds
those x-values to the selected lift, so it rejects the trace. This
explains the failure without reclassifying it as a protocol pass.

The [manifest](evidence/run_36764654520/MANIFEST.json) records seven
GitHub artifact IDs, deterministic compressed archives, every one of the
847 raw file hashes, and the exact run/job conclusions. The
[archive replay](ARCHIVE_REPLAY.json) rehashed all 847 files and compared
the six Mac receipts with the hosted receipts. Five agree in all
cryptographic and decision fields; two last-bit `exp` values differ only
within `1e-12` on n41/L1. The n37 failure reproduces the same assertion.
`archive.py`, `verify_archive.py`, `diagnose_n37_batch.py`, and `analyze.py`
regenerate the receipts and six-row table from the committed artifacts.
For direct raw inspection, make an empty directory and run
`tar -xzf evidence/run_36764654520/artifacts/point-sum-cold-n37_L1024.tar.gz
-C /tmp/point-sum-n37-batch` from this note's directory; the other six
archives use the same form. The manifest must be checked before using an
extracted file as evidence. Second-host receipts retain their exact bytes
in `evidence/run_36764654520/replay/*.json.gz` and can be read with
`gzip -dc`.

The next round should first freeze **new disjoint Q** and an explicit
hosted isolation admission rule that can be met and audited. It should
also represent pinned pair roots as existential x-only witnesses, or
record the exact lifted pair used for a target, then independently verify
both the chosen four points and the pinned root. This correction needs a
new protocol and fresh Q; this run is never selectively repeated.
Given that point-sum is slower than rho in every diagnostic cell and
slower than the S3 batch control at n41, the immediate algorithmic focus
should shift back to useful-base size, relation yield and target PDP
search on n41/n53 batches, with complete cold cost and calibrated
`operations/sqrt(r)` before n83/n131 transfer. The present data do not
support an ECC2K-130 feasibility change.
