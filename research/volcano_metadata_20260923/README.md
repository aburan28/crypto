# Bounded exhaustive validation of volcano metadata

Classification: **accounting / correctness**. This extends PR #614 at parent
`360ebed21e1dc3b637abf2c7c51f36a0cef710c4`. It measures no ECDLP performance;
whole-pipeline operation counts, ratios to rho, and speedups remain null.

The falsification target was zero disagreements with independent finite-field
point enumeration and discriminant predicates, and zero missing kernel edges
between retained vertices inside the expansion depth. Three defects were found
and corrected. Original failed runs are retained.

## Counterexamples and corrections

1. **Incorrect point-count certificate.** For `p=65537, a=1, b=0`, the old
   fast path returned trace 514 instead of 2 (order 65024 instead of 65536).
   `refine_point_order` checked only prime divisors through 97, then treated
   a possibly nonminimal annihilator as the exact point order. It now factors
   the annihilator completely before removing extraneous prime factors. This
   is required for the uniqueness certificate used by the point counter.
   `point-count-failure.log` retains the failure. An additional unit regression
   uses a point of order two and annihilators containing 101, 103 and 127².
2. **Missing edges at the vertex cap.** For `p=5, a=1, b=0`, degree two,
   depth one, vertex cap one, the graph discarded its self-loop. The cap now
   limits discovery of vertices; retained vertices below the depth boundary
   still contribute edges to retained vertices. `initial.log` retains this
   failure. This also preserves backward and parallel edges after the cap.
3. **Broken published evidence link.** The scoreboard's relative research
   link fails after the site builder relocates the page. Replaced it with the
   repository's required absolute GitHub URL. `site-build.log` retains both
   failing assertions; `site-build-fixed.log` records the successful rerun.

## Exhaustive bounds and evidence

The integration test is `tests/volcano_metadata_exhaustive.rs`. Prime fields
are exactly 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43. Every `(a,b)` in each
field is included except singular models; this counts models, not distinct
isomorphism classes. The reference point counter enumerates a histogram of
squares directly, independently of Legendre symbols, BSGS, and point arithmetic.

| Check | Scope | Result |
|---|---|---:|
| Frobenius and conservative endomorphism metadata | Every nonsingular model over the 12 fields | 7,968 passed |
| Fast point-count certificates | Same exhaustive models; `None` is allowed | 5,399 returned, all correct |
| Ordinary/supersingular handling | Within those models | 630 supersingular; 3,820 ordinary orders remain unknown |
| Local position and depth-cap semantics | Degrees 2, 3, 5, 7 excluding the characteristic; caps 0–3 | 31,810 positions passed |
| Discriminant decomposition | Every valid negative order discriminant from -100,000 through -1 | 50,000 passed |
| Rational 2-kernels | Every nonsingular model over the six fields through 19 | 942 models, 876 kernels passed |
| Map evaluation and equal codomain counts | Every affine point for those kernels | 12,936 nonkernel images passed |
| Capped graph edge retention | Every ordinary model through 19; vertex caps 1–4; depth caps 0–2 | 9,984 maps passed |
| Point-count branch boundary | Six fixed coefficient pairs at 16381, 16411, 65537 | 18 counts passed; 17 fast certificates correct |
| Existing CM / volcano / Vélu unit suites | Focused module regressions | 15 / 9 / 5 passed |
| Published site build | Existing offline suite | 39 passed |
| Status snapshot formatting | Existing offline suite | 64 passed |

Discriminant checks reconstruct `D*f²` and independently test squarefreeness
and the fundamental-discriminant congruences. Graph checks retain edge
multiplicity and kernel identity, and compare all eligible edges between
retained vertices. They do not certify completeness beyond the explicit caps.

## Reproduction

```sh
cargo test --test volcano_metadata_exhaustive -- --nocapture
cargo test --lib isogeny::cm::tests
cargo test --lib isogeny::volcano::tests
cargo test --lib isogeny::velu::tests
cargo check --bin crypto
python3 scripts/site/test_build.py
python3 scripts/rho_status/test_rho_status.py
```

Rust 1.98.1 and Cargo 1.98.1 were used. `validation.json` records source/log
hashes and commands. Existing compiler warnings are retained in the raw logs.
The failed GitHub job at the parent commit was the status-pages `test` job;
the site's broken-link assertions were reproduced locally. Local success is
not a claim that a later remote workflow has completed successfully.

## Limits

Exhaustiveness applies only to the stated finite domains. Ordinary
endomorphism claims are tested against the implemented sufficient certificates;
there is no independent general endomorphism-ring oracle. Unknowns are not
converted into guesses. Odd-degree kernel completeness is not established.
These tests neither replay the separate GF(256) Gröbner experiment nor establish
its reported degree of regularity as a curve invariant. They provide no
ECC2K130 relation-yield measurement or evidence of an easier target-size ECDLP.
