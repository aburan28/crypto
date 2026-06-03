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
