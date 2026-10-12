# Round 0006 pre-registration: complete cold batches of 16, continued

Same workload panel as round-0005-batch16: one curve, 16 independent
public-hash targets per cold job, all factor-base, table and precomputation
work charged once to the job, every target solved, every internal check and
the general final scalar check retained, every output verified by the
independent Python checker. Rho solves the same 16 targets through the shipped
per-target signed-Frobenius solver API on the same constructed curve. This
round does not measure the single-target panel, where the last audited result
is 3.04 times rho's instructions and 1.52 times its native time.

Parent: round-0005-batch16 winner `combined_descent` (instruction ratio to rho
0.7195, native 0.7796; parity `true`), restored from the committed evidence
archive and re-verified with its frozen evaluator (1,488 receipts) on this
host before any new work. Source:
`runs/round-0005-batch16/source_candidates/combined_descent/source`. Fresh seed
2026091606; target count 16; `--require-native-progress`; pilot profile (four
development cells, five confirmation cells, 60 fresh confirmation job fixtures,
three repetitions, independent replay). Budget 1,800 paired jobs; seven
challengers and the incumbent schedule 1,680. Single logical CPU 3, 8 GiB
address-space cap, 60-second watchdog per process. Valgrind 3.22.0.

## Host and compiler differ from rounds 0002–0005

This round runs on a different machine (4 vCPU, no SMT siblings exposed)
with rustc 1.94.1, where the archived rounds used rustc 1.98.1 on an AWS host.
Every arm in this round, including the incumbent and rho, is freshly built and
freshly measured under this round's compiler and host, so its ratios are
internally paired. Instruction counts and native times are **not** comparable
across rounds as absolute numbers; only within-round ratios are claimed.

## Why these mechanisms

The round-0005 confirmation profiles were annotated by function (Callgrind,
n23a0, 16 targets). Of the winner's 117M instructions, 58% are shared with the
rho arm (target construction, general final verification, startup, reporting)
and are not touched. Inside the IC-owned 42%:

- `Gf2::inv` (Fermat: `n − 1` squarings and multiplications) is 41% of the
  descent phase, 18% of the setup phase and 57% of the linear-algebra phase.
- `Gf2::new` (the 256-entry reduction tables) is rebuilt once per abscissa
  lift and once per verified relation: 20% of the setup phase.
- Two `FieldStructure` constructions the pair-table path never reads cost
  1.6M each, in the collector and the descent solver.
- `BigUint::modpow` for `λ^k` per relation summand is 13% of descent.
- The Frobenius orbit walks in `finish_factor_base_domain` use big-integer
  squarings and hash keys: about 4.8M of setup.

## Candidates

Seven isolated source changes on the parent, all exact, every check retained;
the last combines five of them. Base support, `m = 3`, batch 4, sparse linear
algebra, the accounting boundary and every verification obligation are fixed.

| id | change | preflight equivalence test |
|---|---|---|
| `it_inv` | Itoh–Tsujii inversion (same power as Fermat, 6 instead of 21 multiplications at degree 23, branch-free) | `inv == inv_fermat` and `a·inv(a) = 1` on every element for `n ≤ 13`, 65,536 elements per larger cell |
| `euclid_inv` | binary extended Euclid over `F_2[z]` | same test |
| `fast_curve_once` | one `FastCurve` per factor-base build and per verified relation batch | lifts and relation verdicts match the per-call construction |
| `lazy_field` | `FieldStructure` built on first use | type-level: same value when read; unread on the pair-table path |
| `lambda_table` | `λ^k mod r` precomputed once per projected orbit map | equals `modpow` for every `k < n` on every cell |
| `fast_orbits` | orbit tables in single-word arithmetic with packed keys | tuple-equal to the general construction on three bases per cell, and identical failure on a non-closed list |
| `combined` | `it_inv` + `fast_curve_once` + `lazy_field` + `lambda_table` + `fast_orbits` | all of the above in one tree |

Expected effect, stated before measuring: the mechanisms together remove at
most roughly 20–25M of the winner's 117M instructions at degree 23, i.e. an
instruction ratio near 0.80 against the incumbent. The 20% promotion gate is
therefore borderline by construction: the shared 58% dilutes every IC-side
gain. A retained incumbent with a measured 10–19% improvement is a plausible
and fully reported outcome, not a failure of the mechanisms.

## Index-calculus admission (new in this round)

Only index-calculus algorithms may compete with rho. Every arm's worker now
reports, for each target, the single relation `[a]G + [b]Q = Σ P_i` over
factor-base points from which the descent derived the logarithm, and the
independent checker verifies that relation in the group and that the returned
scalar is its consequence under the verified column logs. A run whose
logarithms lack that certificate is rejected outright, whatever its cost. The
certificate is emitted by the library descent (`IndividualLogReport.relation`)
and the worker; it is applied to the round's baseline
(`round6-baseline-descent-certificate.patch`) and inherited by every candidate,
so the incumbent is the round-0005 winner plus this reporting. Rounds
0002–0005 ran the same library descent without the per-target certificate. A
first preparation of this round without the certificate was stopped before
any measurement and discarded.

## Reference fairness

`it_inv` and `euclid_inv` change field arithmetic that the shipped rho also
uses. The rho arm is always built from the incumbent's source, so within this
round rho keeps Fermat inversion. Any candidate/rho ratio from an inversion
candidate is therefore provisional: if such a candidate is promoted, the next
round rebuilds rho from the promoted source and that round's rho ratio is the
one to cite. The candidate/incumbent ratio is unaffected.

## Rules

Class: engineering. The floor remains K instructions for K required
independent relation columns; it is weak and cannot establish a non-generic
advance. Within this panel, promotion requires >=20% lower instructions and
native wall, upper paired 95% limits <1 and every cell <=1.10, on confirmation
and replay. Parity requires candidate/rho upper paired 95% limits and every
cell ratio <=1.10 in BOTH metrics on BOTH final stages. No mathematical,
family-wide or cryptographic-size claim follows. Production library defaults
are not changed.
