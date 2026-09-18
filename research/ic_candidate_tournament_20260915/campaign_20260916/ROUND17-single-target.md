# Round 0017 pre-registration: the certificate is the cost — a base named by its orbits

Round 0016 refuted the campaign's strict-win claim on eight cells and handed
forward one hypothesis: that the factor base, sized `6 × degree`, is blind to
subgroup order, and that sizing it by the subgroup would restore the margin at
the two cells where the winner loses to rho outright, `n23a1` and `n29a1`. This
round tested that hypothesis first, at a cost of two minutes, and it is wrong.

## The hypothesis round 0016 handed forward, and why it is wrong

The base is quantised to whole signed Frobenius orbits. `TinyIc::new` draws
representatives eight at a time and only then checks whether the closure has
reached `points`, so every `points` from `4n` to `16n` yields the same eight
orbits, and `20n` to `32n` the same sixteen. Sweeping `points` over the frozen
round-0016 fixtures of the losing cells (four fixtures each, geometric mean of
the per-fixture ratio, `campaign_20260916/round17_base_sweep.py`):

| cell | orbits | `Ir` IC/rho | wall IC/rho |
|---|---:|---:|---:|
| n23a1 | 8 (any `points` ≤ 16n) | 1.2196 | 1.04–1.12 |
| n23a1 | 16 (`points` ≥ 20n) | 1.9636 | 1.27–1.32 |
| n29a1 | 8 | 1.0777 | 0.89–1.02 |
| n29a1 | 16 | 1.8131 | 1.08–1.18 |
| n31a0 | 8 | 0.7606 | 0.89–0.98 |
| n31a0 | 16 | 1.0389 | 1.08–1.14 |

Sixteen orbits is worse at every cell — the pair table grows with the square of
the base and the linear algebra needs twice the relations — and eight is the
smallest size the sampler can produce. The default already sits at the best
accessible size. The instruction ratio is identical across the eight-orbit rows
because the base is identical; the wall varies only by noise, which is itself a
reminder of what the native clock can and cannot see. **Base size is not the
lever.** `n23a1` loses with the same base and nearly the same trial count as
`n23a0` (12 against 11); the cost is somewhere else.

## Where the cost is, cell by cell

Per-phase instructions of the round-0016 receipts for `scan_io` — the round-0015
challenger, which on eight cells already beats rho on instructions in **every
cell** (worst 0.920 at n29a1, 0.881 at n23a1) and panel-wide natively (upper
limit 0.9988), and misses only per-cell native at n29a1 (1.037) and n23a1
(1.021) — medians over 36 receipts per cell, IC-only phases:

| cell | collection | factor base + tables | log certification | IC-only total | rho_solve | `Ir` s/rho |
|---|---:|---:|---:|---:|---:|---:|
| n23a0 | 975,582 | 756,153 | 166,736 | 2,065,820 | 3,003,811 | 0.744 |
| **n23a1** | 1,224,438 | 983,008 | 170,938 | 2,617,925 | 3,132,931 | 0.881 |
| **n29a1** | 77,000 | **878,890** | 208,477 | 1,268,438 | 1,416,625 | 0.920 |
| n31a0 | 597,765 | 903,530 | 237,552 | 1,937,981 | 3,473,230 | 0.604 |

The two losing cells lose for different reasons. At n23a1 the collection is
large — twelve probes at about 100,000 instructions each — and the base build
is a third of that again. At n29a1 the collection is trivial: rho is small
there (`sqrt(42,457) ≈ 206` steps, `rho_solve` 1.42M), and the IC arm's *fixed*
costs — the base build, 69% of everything IC-only, and certification — are
what it cannot amortise. Two independent things nevertheless point at one
object.

First, the instruction-versus-wall gap. At n29a1 `scan_io` is 8% *below* rho on
instructions and 3.7% *above* natively; at n23a1, 12% below and 2.1% above. A
twelve- to fourteen-point decoupling. Round 0016 measured the arm-dependent part
of the process remainder — wall outside the worker's own timer — at +57.7 µs
(n23a1), +77.6 µs (n29a1), +84.4 µs (n31a0) for the IC arm over rho, which is
6–7% of a 1.1–1.3 ms job: the size of the gap.

Second, what that remainder tracks. It is not page faults: rho takes *more*
minor faults than the IC arm (86–95 against 56–87, `getrusage` over seven
runs per cell). Whether it is the evaluator's redirected write of stdout to a
file could not be resolved at this resolution (file-minus-`/dev/null` came out
+115 µs at n29a1 and −20 µs at n31a0, inside the noise). What the remainder
does track, monotonically and at about 5–6 µs per kilobyte, is **output size**:
the IC arm writes 9.7 KB at n23a1, 13.4 KB at n29a1 and 15.0 KB at n31a0
against rho's 0.5 KB, and nearly all of it is the factor base — 368, 464 and
496 points as decimal strings. The same object also costs instructions: the
line-level profile of the round-0016 worker at n29a1 puts the formatting of
those points (`push_u64`, `from_utf8`) at about 250,000 instructions, inside
the phase the dump markers label `log_certification`. That is 13% of the whole
job at n29a1, for writing down 464 points the checker could have regenerated
from eight.

## What changes, and for whom

**`orbits`** — the round-0016 `scan_io` source plus
`campaign_20260916/round17-orbits.patch`: the certificate names the factor base
by its **orbit representatives** instead of listing every point. One point per
signed Frobenius orbit (eight at every cell but n13a0, which closes to seven),
in the worker's orbit order, under a new key `factor_base_orbits`. Nothing else
in the certificate changes: the relations, the column logs and the descent
relations index the same flat base at the same positions.

**The checker** — `campaign_20260916/round17-oracle-orbits.patch` to
`oracle.py`, additive: a report carrying `factor_base_orbits` has its flat
base regenerated by the oracle's *own* `frob` and `neg`, as `+R, −R, +Frob(R),
−Frob(R), …` for each representative `R`, with the orbit required to close
under Frobenius; a report carrying `factor_base` verifies exactly as before;
a report carrying both is refused. Every downstream check — relation
membership in the group, the scalar-field rows, rank, descent — runs on the
regenerated list, and `factor_base_sha256` is over the same points. The format
was admitted only after two checks on round-0016's frozen reports: expanding
the representatives of 24 `scan_io` reports across all eight cells reproduces
their listed bases byte for byte, and verifying 16 reports in both formats
gives identical base hashes, relation counts and ranks, while a report whose
first representative has its coordinates swapped is refused ("point does not
lift to curve").

This is a protocol amendment, and it is stated as one: the checker gains a
second accepted certificate format. It is smaller, not weaker — the checker
recomputes every point it uses from the representatives with its own
arithmetic, and a wrong representative fails every relation that touches its
orbit. The size of the certificate shrinks by `2n`: 13,422 to 206 bytes at
n29a1.

What this round deliberately does **not** change, and why:

- The serialisation itself. `scan_io` already formats digits two at a time
  into a byte buffer (round 0015's `fastio`); what remains is volume, which
  this round removes.
- The second cofactor multiplication per orbit. `TinyIc::new` projects each
  sampled point into the subgroup with `[h]` and then forms the column as `[h]`
  of the representative again — at n29a1, where `h ≈ 12,644`, sixteen 14-bit
  scalar multiplications for about 275,000 instructions, 30% of the base
  phase. Dropping it would change the column definition, which the checker
  hard-codes (`[h]P` locates a base point's column; rows check `h·a`), so it
  is a second checker amendment. It is left for a round of its own, with its
  own equivalence check, rather than folded into this one.
- The collection at n23a1. Twelve probes at ~100,000 instructions is the
  scan round 0015 already attacked; nothing in this round touches it, and
  the pre-registration below says what that means for n23a1.

## Objective and gates

Objective `rho`, unchanged. Promotion over the incumbent requires the
instruction ratio at or below 0.98 with a 95% upper limit below 1, a native
95% upper limit below 1, and every cell within 1.10, on confirmation and on
replay. `beats_rho_strict` requires the winner over rho below 1 in both metrics
in every cell, on both stages. No gate, margin, bootstrap, repetition count,
confidence level or panel changes in this round. The one protocol change is the
checker's second accepted certificate format, stated above.

## Parent, arms, seed, budget

Parent: round 0016 (`runs/round-0016`). Three arms: `incumbent`, the retained
winner exactly as sealed there; `scan_io`, the round-0015 challenger exactly as
sealed there, carried as a **control** so that any difference between it and
the new arm is the certificate format alone; and `orbits`, `scan_io` plus
`round17-orbits.patch`. Panel: the round-0016 eight cells — `13a0, 17a1, 19a0,
23a0, 23a1, 31a0` with `19a1, 29a1` held out to confirmation and replay. Seed
2026091717, pilot profile, one target per job, CPU 3, 60 s timeout,
`--max-processes 3000` (2,268 trials: 36 aa, 72 smoke, 216 development, 216
selection, 864 confirmation, 864 replay).

## Boundary, floor, class, honesty

Unchanged from round 0016. The unit is Valgrind 3.22 `Ir` over the complete
worker process; native timings are paired cold-process wall under the same
blocking-reap protocol and carry their own scope. The class is accounting, not
a complexity claim, and nothing here is a statement about any curve outside
the eight cells named above. A smaller certificate is not a weaker one here —
the checker regenerates every point it verifies against — but it is a
different one, and the amendment that admits it is recorded as such.

## What this round expects, stated before any tournament stage

Development measurement of `orbits` against `scan_io` on four round-0016
confirmation fixtures per cell (`orbits` verified to the identical certificate
on all 96 fixtures first; instructions by callgrind, wall the median of five
interleaved cold runs; `campaign_20260916/round17_measure.py`):

| cell | `Ir` orbits/scan_io | wall orbits/scan_io |
|---|---:|---:|
| n13a0 | 0.938 | 0.979 |
| n17a1 | 0.924 | 0.921 |
| n19a0 | 0.933 | 1.065 |
| n19a1 | 0.943 | 0.987 |
| n23a0 | 0.959 | 0.974 |
| **n23a1** | 0.966 | **1.005** |
| **n29a1** | 0.912 | **0.889** |
| n31a0 | 0.934 | 0.939 |

Predictions, in the order they will be checked:

1. **`orbits` is byte-faithful.** On every confirmation and replay fixture its
   relations, column logs and solutions are identical to `scan_io`'s and its
   base hash is identical; only the certificate's length differs. Any
   difference is a bug in the patch and ends the round.
2. **`orbits` over `scan_io` is a larger native gain than instruction gain at
   n29a1, and no native gain at n23a1.** Development says −9% instructions
   against −11% wall at n29a1, and −3% against 0% at n23a1. If the remainder
   round 0016 found tracks output size, confirmation will show the same shape:
   `orbits`/`scan_io` native below its instruction ratio at n29a1, and near one
   at n23a1. If instead the native ratio simply follows the instruction ratio
   everywhere, the output-size explanation is wrong.
3. **n29a1 flips.** `scan_io` was 1.037 over rho natively there in round 0016;
   at 0.889 of `scan_io`, `orbits` is expected near 0.92 — below one with
   margin against the ±5% per-cell spread the A/A control shows.
4. **n23a1 does not flip, and `beats_rho_strict` therefore still fails.**
   `scan_io` was 1.021 over rho natively at n23a1 and `orbits` gains nothing
   there natively; `orbits_rows` is predicted to cut n23a1's instructions by
   about 2.7%, which is not enough on its own. This round expects the strict
   gate to fail at exactly one cell, n23a1, on native wall only, with every
   cell below one on instructions. That is stated now so that a strict pass
   would be a genuine surprise and a fail at any *other* cell a real anomaly.
5. **`orbits_rows` stores two rows at n23a1 and the same rows as `orbits`
   everywhere else**, by the weighted rule's own arithmetic; its instruction
   ratio to `orbits` is below one at n23a1 and one elsewhere.
6. **Promotion.** `orbits` is expected to pass the no-regression gate over the
   incumbent on both stages: its instruction ratio is about 0.83 and its
   native upper limit, which `scan_io` missed by 0.0034 on replay, falls with
   the certificate. If it does not pass, the reason will be a native upper
   limit again, and the pre-registration of round 0016 already says why that
   gate cannot resolve small effects.

Falsification, stated before the run: if prediction 1 fails, the round is
void. If the strict gate fails at any cell other than n23a1, or on instructions
anywhere, the output-size account of the remainder is not the whole story and
the decision record says so. If it passes, the campaign's rho claim is restated
at eight cells — and that restatement, unlike round 0013's, will have survived
the panel that refuted its predecessor.

What comes next regardless: n23a1 loses on the collection — twelve probes of
about 100,000 instructions each — and no format change reaches it. Round 0018,
if this round lands where it is predicted to, is about that collection.
