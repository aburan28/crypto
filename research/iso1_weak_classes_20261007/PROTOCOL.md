# ISO-1 trace census protocol (documented 2026-10-07)

This file records the implemented method and stop rules. It was written
after the p = 37 run and after the p = 41 run had started; it is **not**
a claim of prospective registration for those runs. The original user
request is p = 37 through about 200, q = p², every trace, with a
necessary-and-sufficient class criterion.

## Hypotheses and success conditions

1. The correction to the second F_(p²) square-class branch changes
   historical weak-class counts. Success: use a genuine F_(p²)
   nonsquare and enumerate 2q²+2q normalized representatives for each
   prime. The full and twist-derived methods must agree at a small p.
2. Frobenius conductor depth at 2 predicts weak-class membership.
   Success as a **necessary condition**: no observed weak trace with
   v₂(f_pi)=1. Success as an **exact formula**: no false positives or
   false negatives on every trace at every requested p, plus a proof
   covering arbitrary p. A held-out error leaves the exact-formula
   hypothesis rejected.
3. The sampled reach fraction remains in the historical [0.50,0.68]
   band at p=37,41,43. Success for a sampled prime: 4,000 random
   full-2 curves, fixed `StdRng` seed `1 ^ 0xE8AC7`, with a Wilson
   interval and the corrected weak trace set. This is a sample of
   curves, not a complete count of all curve isomorphism classes.

## Frozen construction and trace record

Use `examples/iso1_class_census.rs` and `Fq3::new(p)` from the checked
source tree. For p≠3, F_(p²)=F_p[u]/(u²−w), with first F_p nonsquare
w, and F_(p⁶)=F_(p²)[theta]/(theta³−s), with first noncube s. For
every a0∈F_(p²), enumerate alpha=a0+a1 theta+a2 theta² with
`a1∈{1,w2}` and `a2∈F_(p²)`, plus `a1=0` and
`a2∈{1,w2}`; `w2` is the first F_(p²) nonsquare. That gives
2q²+2q normalized representatives. In `--derive-twists` mode,
evaluate one square-class branch and derive the other by trace negation.

Write one CSV row for every Hasse-interval `t≡2 (mod 4)`, including
zero-witness and nonordinary rows. Record representative count,
fundamental discriminant, Frobenius-order conductor depth, 2-splitting,
v₂((t/2)²−p⁶), and maximal-order class-number parity. Treat
`curve_order` as a probabilistic point counter; a completed enumeration
is not an independent certificate for every assigned trace. Keep
uncertain labels explicit and retain failed audit evidence.

## Limits and accounting

Use `RAYON_NUM_THREADS=4` for large twist-derived runs. Save the
command, exit code, binary hash, wall seconds, and charged F_p
multiplications. This is a diagnostic count, not a paired IC/rho
benchmark. Stop a run on an actual crash, disk exhaustion, or a
documented explicit resource limit; preserve the failure. Do not infer
that the exact formula holds because a sampled or bounded run finds no
counterexample. A proposed formula cannot replace the cover reach
census or seed sieve until it is both sufficient and proved.

The current direct algorithm grows roughly as p^4 representatives
times p^(3/2) point-count work. The measured p=37 run and the
extrapolated p=199 cost are in [REPORT.md](REPORT.md). This is a
computational limit of this implementation, not a mathematical
impossibility statement.

## 2026-10-08 amendment: exact orbit quotient

After the direct p=37,41,43 runs, the norm-one proof gave an exact
sixfold orbit reduction. The projective normalized alpha values map
bijectively to nonidentity elements of the norm-one torus by
`lambda=alpha^q/alpha`. The transformations `lambda -> lambda^q` and
`lambda -> lambda^-1` preserve the symmetric trace pair. They have
six-element orbits except the two nonidentity cube roots of unity,
which form one two-element orbit. `--derive-twists --orbit-quotient`
selects the least encoded lambda in each orbit and weights its trace
by six or two. The exact expected point-count call count is
`(q²+q+4)/6`; the existing twist-derived count is `q²+q`.

Validation before using this mode for further primes: compare its
entire CSV byte for byte against a completed twist-derived run at a
small p; check the representative and point-count call totals; retain
any differences from historical full-mode data, including zero-label
changes. This mode inherits the randomized point counter's uncertainty.
Its operation-count reduction is exact; wall-time performance claims
would require the repository's CPU isolation gate.
