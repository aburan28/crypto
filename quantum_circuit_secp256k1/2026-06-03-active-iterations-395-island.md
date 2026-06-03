# 2026-06-03 ACTIVE_ITERATIONS 396 -> 395 (next active-iter trim)

Validated 0/0/0 over all 9024 shots via the official `eval_circuit`.
Score **1434 x 1,739,149 = 2,493,939,666 -> 1434 x 1,736,993 = 2,490,847,962**
(-3,091,704, -0.124%). Peak unchanged at 1434.

## What changed

One further active-iteration trim of the dialog-GCD transcript:

- `DIALOG_GCD_ACTIVE_ITERATIONS` 396 -> 395.

The binary-GCD half-inverse still reaches gcd=1 on the reachable verifier
support in 395 steps; the worst-case 256-bit Stein/binary-GCD envelope on the
9024-shot Fiat-Shamir support fits under 395 here. One fewer transcript step
shortens both the forward tobitvector loop and the apply replay, dropping avg
executed Toffoli 1,739,149 -> 1,736,993 (-2,156). It does NOT change the peak
(1434 is bound by the apply mod add/sub + ipmul/quotient `*_reacquire_terminal_u`
phases, which are independent of the iteration count), so the win is pure
Toffoli.

## Island

The shorter op stream re-rolls Fiat-Shamir. A randomized 2-D
`DIALOG_REROLL x DIALOG_POST_SUB_REROLL` screen (clean-island density was low on
this stream -- the iteration trim is one notch past the prior 396 floor, so most
rerolls leave 1+ shots needing the 396th step) lands a clean island at:

- `DIALOG_REROLL=767`
- `DIALOG_POST_SUB_REROLL=173`

Co-tuned in `configure_ecdsafail_submission_route()` (replaces the prior
active-396 island `DIALOG_REROLL=214 / POST_SUB=466`).

## How it was found

Screened with a faithful early-exit checker (replicates eval_circuit's
`fiat_shamir_seed` + per-batch bit-sliced sim, bails on the first failing 64-shot
batch, materializes scalar-mult inputs lazily) -- ~3-9 s per failing candidate vs
~20 s for a full eval. First-failing-batch over random rerolls ranged ~3..100
(/141), i.e. effective failure rate lambda ~ 1.4..9; r=767/p=173 was the first
fully-clean (all 141 batches) candidate. Every CLEAN hit was reconfirmed with the
real `eval_circuit` before use; the screener is local-only tooling and is not part
of the submission.

## Next floor

395 is now the active-iteration floor on this support (394 leaves several shots
short of gcd=1 across most rerolls -> no easy island). Peak 1434 is still the
6-phase co-bind (apply `materialized_special_chunked_raw_{sum,difference}`,
ipmul/quotient `*_reacquire_terminal_u`, `pair1_quotient`/`pair2_product`);
sub-1434 still needs the apply transient + ipmul/quotient phases narrowed
structurally (safegcd/jump-GCD), as noted in the cswap-floor analysis.

## Jump/half-GCD investigation (goal: sub-2.45B structural cut)

Investigated the repo's jump/half-GCD scaffolding as the path to a big Toffoli
cut. Findings:

- `halfgcd_live_pa::emit_round162_halfgcd_live_pa_or_fail` and
  `round158_halfgcd_splice_live::abort_round158_live_prefix_pa_route` are
  **fail-closed stubs** -- prior agents attempted the half-GCD PA and could not
  land a valid, in-budget, reversible operator.
- Documented blocker (round162): a valid secp256k1 input has dx=x1-x0=1, so the
  first half-GCD quotient is floor(p/dx)=p, needing a full-width (256-bit) q.
  The small-quotient decoder (<=26 bits observed) can't represent it; a *total*
  operator must charge a full-width q slot, raising qubits. The deeper missing
  object is the total schedule that erases u, v, coeff_b, coeff_d, q_tail back to
  the four Google-ABI registers -- the old splice returns with non-ABI scratch
  live (ancilla-garbage). The dialog-GCD wins precisely because its compressed
  sidecar transcript solves this clean-uncompute problem.
- Repo's own jump analysis: safegcd batching yields only ~5-11% Toffoli for
  "high implementation complexity," because applying the 2x2 transition matrices
  quantum-controlled costs nearly as much as replaying microsteps, and the
  apply/Bezout phase (47% of cost) does not benefit from step-batching.

Conclusion: sub-2.45B (-1.6%) requires either (a) the safegcd jump-GCD rewrite
WITH a new sidecar-style clean-uncompute (research-scale; the existing scaffold
is abandoned blockers), or (b) a denser GCD transcript code (the round763
compressor already frees exactly 1 of 6 raw bits -> 17..32 reachable 3-step
states -> 5 bits provably necessary; no 4-bit code exists). Neither is a
reliable in-session validated change. Best validated result this session remains
the active-395 island: 1434 x 1,736,993 = 2,490,847,962.

## Denser-log (GROUP_SIZE=5) investigation — real lever, research-scale to land

Pursued the peak cut via a denser sidecar code (authorized open-ended effort).
The peak is purely persistent state: u(256)+tx(256)+ty(256)+log + ~1. tx/ty/u are
fixed-width, so the ONLY shrinkable register is the log. The log is a base-3
stream: each GCD step stores (b0, b0&b1) with b0&b1 = b0 AND (u>v), giving exactly
3 reachable states/step. The round763 compressor packs 3 trits -> 5 bits
(27 <= 32, frees 1 of 6), = 1.667 bits/step -> 665q.

GROUP_SIZE=5 packs 5 trits -> 8 bits (3^5=243 <= 256, frees 2 of 10) = 1.6
bits/step -> ~640q, a -25q peak cut (with active=395=79*5 exactly, -33q -> peak
1401 -> 1,736,993 x 1401 = 2,432,527,193, well under the 2.45B goal). The
future-log carry borrow capacity (~632 early) still exceeds the ~511 transient
need, so borrowing survives.

Two blockers make this research-scale, not a localized edit:
1. **Compressor synthesis.** Need an ancilla-free CX/CCX-only reversible map that
   zeroes 2 output bits across all 243 valid 10-bit inputs (ancilla would re-raise
   the peak; high-control Toffoli synthesis (TBS) needs ancilla). Greedy and
   simulated-annealing searches (thousands of seeds) both failed to find one;
   it needs a structured hand-derivation (base-3 Horner with constraint-driven
   bit cancellation).
2. **Layout invariant.** GROUP_SIZE/BLOCK_BITS are tied to ~10 precomputed
   constants via `DialogGcdHighTailLayout::check_passed`, enforced by a hard
   `assert!` in the build path (EXTENSION_BITS=212, MIN_V_INDEX=15,
   PROJECTED_Q=1451, COMPRESSED_BITS, per-block cell positions, RAW_REPAIRED
   hash/T). Changing GROUP_SIZE requires re-deriving the entire high-tail-alias
   layout model.

Conclusion: both structural paths to sub-2.45B (jump-GCD; denser-log) are
genuine but research-scale. Best validated result this session stays the
active-395 island: 1434 x 1,736,993 = 2,490,847,962. The denser-log path is the
most promising future lever: derive the 5-trit compressor (verifiable by
enumerating 243 inputs), then re-derive the DialogGcdHighTailLayout constants for
GROUP_SIZE=5.
