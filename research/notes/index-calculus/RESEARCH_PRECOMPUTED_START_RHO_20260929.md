# Precomputed-start strong rho: does the one favourable compact-orbit cell survive?

Ledger rank 1 (`docs/ic/boundary_targets.json`, `agent_priorities[0]`) asks for a
precomputed-start signed-Frobenius rho before any further `vs_rho` claim. The
strong-rho sweep (`RESEARCH_STRONG_RHO_SWEEP_PROTOCOL_20260929.md`) already found
IC/rho ≥ 1.5 in retired instructions at every cell with n ≥ 41, rising with n. The
only cell below 0.8 is n = 37, L = 1,024, K = 42 at **0.794**, which that note
attributes to fixed costs, among them rho's per-target start scalar
multiplications. This note tests exactly that cell (and n = 41 as a control)
against a rho whose starts are precomputed.

## What was already seen before writing this (disclosure)

The x86-64 sweep table (IC/rho R3 = 0.794, 1.520, 1.923, 2.297 at n = 37, 41, 53,
61; L = 1,024) and the ladder note. **Not seen:** any rung-4 number, or any
number for any arm on the host used here.

## Rung 4 (new file `examples/koblitz_rho_batch_ks_strong_ps.rs`)

A copy of `examples/koblitz_rho_batch_ks_strong.rs` (unchanged, so the ladder's
hashed binary is untouched) with one added rung. R4 = R3 (library `Gf2`,
least-rotation canonicalization, 32-lane lockstep batch inversion) plus:

- the walk stride `[s]G` and a base start `[c]G` are computed **once** in setup
  (two scalar multiplications, charged), from a deterministic seed separate from
  the jump table;
- for each target the first start is `[c]G + Q` (one charged addition) and later
  starts advance by the stride (one charged addition each), exactly as in R3;
- the per-target stride and cursor scalar multiplications of R3 are removed.

Unchanged and still charged in both arms: target generation `Q = [d]G` from the
public fixture (IC also generates Q from `scalars.txt`), collision verification
`[candidate]G = Q`, jump table, distinguished-point table. Starts are
`[c + k·s]G + Q_t`, so walks of different targets never share a start point.

## Cells (frozen; identical to the sweep)

| n | L | K (IC) | corpus | batch_seed | dp_bits | rank_seed |
|--:|--:|--:|:--|--:|--:|--:|
| 37 | 1,024 | 42 | `n37-strong-sweep-L1024-v1` | 531310 | 4 | 7 |
| 41 | 1,024 | 255 | `n41-strong-sweep-L1024-v1` | 531310 | 4 | 7 |

Commands as in `strong_rho_sweep_20260929_run/sweep_cell.py`: rho
`KIC_RHO_RUNG={3,4} KIC_RHO_LANES=32 KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=<corpus>
koblitz_rho_batch_ks_strong_ps <n> 0 signed_frobenius 1024 531310`; IC
`koblitz_orbit_dlp_fast construct:<n>:0:<K> scalars.txt 7 <out>`, with
`scalars.txt` the rho corpus' planted scalars.

## Host and metric (fixed now)

One Apple-silicon macOS host (arm64), release build, single-threaded processes,
run only while the 1-minute load average is below the core count's machine limit
(14). Metric: **"instructions retired"** reported by `/usr/bin/time -l` for the
whole process, five repeats per arm per cell, alternating arms, median reported.
Valgrind does not run on macOS arm64, so these numbers are **not comparable** to
the x86-64 callgrind `Ir` of the sweep; only ratios within this host are used.
Wall, user+sys CPU and peak RSS are recorded as secondary.

## Gates

- **G1.** Every arm: all targets recovered and verified; rho equal to the planted
  scalars; IC zero failures and full rank.
- **G3.** R4 total walk steps within `[0.85, 1.25]×` R3's on the same cell.
- **G4.** Unit test: R4 recovers every target end to end at n ∈ {17, 19, 23}.

## Decision rule (fixed now)

Let `ρ₃ = IC/R3`, `ρ₄ = IC/R4` (medians of instructions retired, this host), and
`ρ* = IC / min(R3, R4)` over rungs passing their gates.

- Host transfer: `ρ₃` at n = 37 is reported against the x86 0.794 without a verdict
  (different ISA and unit).
- **The favourable cell survives precomputed starts** on this host if `ρ* < 0.8`
  at n = 37; **it does not survive** if `ρ* ≥ 1.0`; `[0.8, 1.0)` is reported as
  inconclusive.
- n = 41 is a control: `ρ*` there is expected ≥ 1 from the sweep; a value < 0.8
  would contradict the sweep and is reported as such, not explained away.

**Inadmissible, stated in advance:** changing L, K, dp_bits, lanes, seeds or
corpora; dropping target generation or verification from either arm; scoring on
wall clock; reporting a rung that failed a gate.

## Scope

Two cells, one arm64 host, IC as merged (no comparable optimisation), R4 not the
best possible rho. Nothing here speaks to n ≥ 53 beyond the sweep, to other L, or
to n = 83 / ECC2K-130 transfer. Public synthetic known-answer fixtures only; no
key recovery.
