# Archived multi-target IC panels against rho

> **Measurement status:** the tables below are explicitly multi-target batch
> experiments. They do not answer the primary one-target online question and
> must not be quoted as single-target speedups. The current primary paired
> one-point run is documented at
> [`experiments/koblitz-single-target-compliance-20260924/README.md`](../../../experiments/koblitz-single-target-compliance-20260924/README.md):
> one IC target and one rho walk per invocation, with target-independent IC
> preparation excluded from online time. The batch panels remain secondary
> results for the separate shared-log/amortization question.

On matched public-synthetic n=53 targets, the compact-orbit shared-factor-log DLP
beats frozen Kuhn–Struik batched rho in all 15 paired n=53 blocks tested
(L = 32 to 16,384), and in all 3 blocks at n=41 (L=1,024). It is 1.3–1.8x faster at L ≤ 1,024 with the retained K=220 base,
and about 4x faster at L = 4,096 and 16,384 with a larger base chosen by a
pre-declared rule. It has no pair table and no edge selectors. It uses 0.36–5.7 GB,
against batched rho's 0.06–0.78 GB. Every target's relation and recovered discrete
log was replayed in independent pure-Python arithmetic with zero failures.

This is a finite constant-factor result. At their optimal parameters both arms
cost about √(L·r/n), so it is not an asymptotic sub-rho claim.

## What changed

The previous compact extractor (`koblitz_s5_sat_instance`, SAT-confirmed) needed a
median 667 ms per n=53 relation and rebuilt its index for every target.
`examples/koblitz_orbit_dlp_fast.rs` changes four things:

1. **One index per run.** The Frobenius-quotiented S3 regular-root table is built
   once and shared by the rank stage and every target.
2. **One probe per partner.** In a normal basis, Frobenius is a bit rotation. Each
   root is stored under its smallest rotation, so a demanded partner x needs one
   table probe instead of up to n squarings and n SipHash lookups. Squaring,
   half-trace and trace are byte-table linear maps. Inversion uses Itoh–Tsujii
   with the `a^(2^k)` steps done as rotations.
3. **Rank-guided relations.** Decomposing `[a]G − R_j` for a pivotless column j
   makes every row raise the rank, so the rank stage needs exactly K relations
   (220 at n=53). Random relations needed 413–476.
4. **Base size K as a parameter.** The base can be built in-process by a public
   x-scan (charged to the run). Relation cost falls as 1/K² while index cost grows
   as K², so the best K grows with L.

Each query scans from a target-dependent state offset. A fixed start sent every
first hit to the same early columns.

## 1. Growing n against independent rho (retained bases, 32 targets per block)

All stages run and are timed in one process: base load, index, rank stage,
linear algebra, targets and verification. The 12 blocks are the 2026-09-12
panel's fixtures and bases.

| n | log₂ r | Compact total (median) | Explicit-table IC | Independent rho | Compact / rho | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|
| 37 | 27.78 | 25 ms | 38 ms | 185 ms | 0.134 | 8 MB |
| 41 | 39.00 | 917 ms | 2,135 ms | 5,193 ms | 0.175 | 55 MB |
| 53 | 44.26 | 8,816 ms | 28,262 ms | 37,187 ms | 0.255 | 359 MB |

Slope on log₂ r: compact 0.504, explicit table 0.571, independent rho 0.459. At
fixed K per n the compact arm still grows faster than rho. All 384 targets were
recovered and replayed.

## 2. Paired against Kuhn–Struik batched rho at n=53

Protocol: `paired_vs_batched_rho_protocol.json`. Batched rho is the frozen
`autolab_batched_rho_n53_20260922/frozen/v2` binary with d=4, as in that lab's
panel. It generates each corpus, and the compact arm receives the same published
scalars. Both arms are single-threaded. There are 3 blocks per L, and block parity
alternates which arm runs first. The decision rule is the batched-rho lab's: an
arm wins at L only if all 3 block wall ratios favour it. The host was shared, with
load averages 27–40.

| L | K | Compact wall | Rho wall | Wall ratio (range) | User-CPU ratio | Peak RSS compact / rho | All 3 blocks |
|---:|---:|---:|---:|---:|---:|---:|---|
| 32 | 220 | 9.3 s | 12.8 s | 0.742 (0.729–0.746) | 0.716 | 359 / 58 MB | compact |
| 1,024 | 220 | 35.7 s | 63.9 s | 0.559 (0.542–0.586) | 0.553 | 359 / 212 MB | compact |
| 4,096 | 220 | 112.6 s | 122.3 s | 0.912 (0.908–0.934) | 0.924 | 359 / 370 MB | compact |
| 4,096 | 660 † | 30.2 s | 116.4 s | 0.243 (0.239–0.259) | 0.233 | 2.7 / 0.37 GB | compact |
| 16,384 | 880 † | 63.7 s | 243.6 s | 0.265 (0.262–0.306) | 0.241 | 5.7 / 0.78 GB | compact |

† K was set before the run by fitting T(K) = aK² + (K+L)·b/K² to a K = 220/440/660
sweep at L=1,024 on the panel-1 corpus. The constructed base hash is identical
across blocks at each K.

With K=220 fixed, the gap closes as L grows (0.56 at L=1,024, 0.91 at L=4,096).
Compact cost per target is flat, while batched rho's cost per target falls with L.
Raising K with L restores the gap.

Every fixture's published point and recovered scalar agree between arms, across
all 15,456 (panel 1) and 61,440 (panel 2) targets.

## 3. Growing n against batched rho at fixed L = 1,024

Script: `growing_n_vs_batched_rho.sh`. For each n, K was chosen by lowest wall on a
disjoint tuning corpus (`n<n>-ks-growing-tune-1024-v1`) before the evaluation
corpus (`n<n>-ks-growing-1024-v1`) was run. There are 3 paired blocks per n,
alternating arm order. Load averages were 28–33.

| n | log₂ r | K | Compact in-process | Batched rho in-process | Wall ratio (range) | Peak RSS compact / rho | All 3 blocks |
|---:|---:|---:|---:|---:|---:|---:|---|
| 37 | 27.78 | 42 | 360 ms | 401 ms | 0.950 (0.848–1.079) | 14 / 7 MB | unresolved |
| 41 | 39.00 | 255 | 2,717 ms | 8,606 ms | 0.339 (0.276–0.358) | 350 / 32 MB | compact |
| 53 | 44.26 | 440 | 14,926 ms | 59,513 ms | 0.246 (0.242–0.264) | 1.35 / 0.21 GB | compact |

At n=37 both arms take about 0.4 s and fixed process setup dominates. The in-process
slope on log₂ r is 0.467 for compact against 0.531 for batched rho over n = 41→53,
and 0.316 against 0.431 over all three sizes. These fits use 2–3 points with tuned
K and have no interval. The √(L·r/n) optimum is the same for both arms, so read the
difference as finite-size, not as a growth-rate claim.

The next usable Koblitz rung above n=53 is a=0 n=61 (r = 162888033982417, 48 bits,
cofactor 14,156). a=0 n=59's remaining 57-bit factor is composite. The n=61-admitting
rho binary (`examples/koblitz_rho_batch_ks_v2_n61.rs`) differs from frozen v2 by one
admitted-n line and is bit-identical at n=53 on a matched corpus.

**Do not run 32-target panels.** The next measurement is n=61 vs batched rho at
**L ≥ 1,024** (K tuned on a disjoint 1,024-target corpus; then L=4,096 if RSS
allows). Script: `growing_n_n61_L1024.sh`.

## Independent replay

`independent_replay.py` uses its own pure-Python GF(2^n) and Koblitz arithmetic,
with no code shared with the producer. For each record it checks:

- the target is [scalar]G and has order r,
- every x-code is in the base,
- some sign choice of the four base points sums to the target,
- the pinned S3 intermediates match,
- every used base point has order r,
- the recovered log d satisfies [d]G = target and equals the published scalar.

| Set | Records | Pass | Fail |
|---|---:|---:|---:|
| Matched slope panel (SAT-confirmed, n=37/41/53) | 384 | 384 | 0 |
| Fast DLP vs independent-rho fixtures | 384 | 384 | 0 |
| Paired panel 1 | 15,456 | 15,456 | 0 |
| Paired panel 2 | 61,440 | 61,440 | 0 |
| Growing n vs batched rho | 9,216 | 9,216 | 0 |

## Limits

- The host was shared, with load averages of 27–50 from other sessions. Load was
  recorded per run, and user-CPU ratios are reported alongside wall.
- There is one implementation of each arm, on one host, at n=53 for the batched
  comparison.
- Batched rho is the frozen v2 producer with d=4. Its memory is 7x lower at
  L=16,384. A Bernstein–Lange precomputed table, or a better-tuned batched rho, was
  not compared here.
- At optimal K both arms scale like √(L·r/n), so this is a constant-factor
  comparison. At fixed K, the compact arm's growing-n slope is above rho's.
- Replays were run by the same agent on the same host. Rerunning the producer on
  an independent host is still pending.
- Public synthetic fixtures only. No key recovery, external points, deployed
  curves or asymptotic claim.
