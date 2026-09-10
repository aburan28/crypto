# Stage 16: source-equivalent compact GGMP Magma evidence

The verbose Stage 13 Magma exports repeat variable expressions and put every `n=31` GGMP file above the official calculator's 50,000-byte cap. Stage 16 changes representation only: it parses each WDSat ANF row into Boolean monomial masks, independently parses each archived Magma `F` polynomial into the same masks, requires exact ordered monomial-list equality, and renders Magma's compact `BooleanPolynomial(R, Q)` form. A retained three-variable Magma sanity run confirms the integer-mask orientation by comparing verbose and compact polynomials inside Magma.

This equivalence check passes for all twenty frozen Stage 13 tasks. The compact renderer retains the same Boolean ring rank and `grevlex` order, equation sequence, `SetGPU(false)`, deterministic seed, direct Faugère F4, sparse mode, and `Nthreads := 1`. Four compact inputs fit the observed service cap: the two Stage 15 standard cells plus GGMP seeds `2026091303` and `2026091305`.

| Seed | Cell | Compact bytes / SHA-256 | F4 result | Step degrees | Basis size | Internal CPU-s | Internal wall-s | Service time-s | Service memory field |
|:--|:--|:--|:--|:--|--:|--:|--:|--:|:--|
| 2026091303 | n31 GGMP | 47,680 / `99a216bd87d9d77007120ee9a0d907a25d0186ae103c607f0f88a1de5965a164` | proper ideal / SAT certificate | `[2,2,3,3,4,3,4]` | 87 | 13.760 | 13.770 | 13.779 | `267.94MB` |
| 2026091305 | n31 GGMP | 48,284 / `73237e89419aa7dc89e06bab33cc6093af0f381b494549eeee40e29402f6e646` | proper ideal / SAT certificate | `[2,2,3,3,4,3,4]` | 87 | 13.820 | 13.820 | 13.830 | `267.88MB` |

The raw XML contains exactly the seven expected F4 terminal markers and no service warning. Both calculator-reported terminal identities agree on algorithm, status, degree sequence, basis size and controls. The basis polynomials were not exported and no model was extracted, so both remain `sat_basis_certificate_unverified_model`.

The service time and memory fields are retained as opaque service metadata, not process-scoped CPU or peak RSS. The internal timers surround only the direct F4 call. These values must not be ranked directly against the Stage 13 process envelopes.

This extends supporting Magma coverage to two source-equivalent standard and two source-equivalent GGMP `n=31` systems. It does not execute `n=41` or `n=59`, supply five replicates per cell, validate Magma-linked rational point witnesses, close resource accounting, or change the non-SOTA conclusion.

The compact constructor and Boolean-ring semantics are documented in the official [Magma Gröbner-basis handbook](https://magma.maths.usyd.edu.au/magma/handbook/text/1314). The retained verifier and summary bind the compact inputs to the exact Stage 13 ANF and verbose Magma exports.
