# Stage 17: bounded n=41 Magma F4 evidence

A named-generator representation reduces the source-equivalent Magma input size without changing the Boolean polynomials. Each source variable is assigned one distinct Boolean-ring generator; every term is rendered from the same ordered monomial masks independently recovered from the frozen ANF and verbose Magma exports. The Stage 16 all-twenty equivalence and Stage 13 custody checks run before this receipt is admitted.

The named representation places all ten `n=31` and all five `n=41` tasks below the official calculator's observed 50,000-byte cap. All five `n=59` tasks remain 182,417 to 195,340 bytes. To obtain a bounded literature-scale check without running a bulk public-service campaign, Stage 17 submitted only seed `2026091301` of the `n=41,ell=5,m=3` standard cell.

| Seed | Cell | Named input bytes / SHA-256 | F4 result | Step degrees | Basis size | Internal CPU-s | Internal wall-s | Service time-s | Service memory field |
|:--|:--|:--|:--|:--|--:|--:|--:|--:|:--|
| 2026091301 | n41 standard | 22,517 / `33bb9cdc0d54247009978d889df1b1d69118e73414011951615bb090b9c973a6` | proper ideal / SAT certificate | `[2,3,3,3,3]` | 43 | 0.740 | 0.730 | 0.780 | `32.09MB` |

The raw Magma V2.29-10 XML contains the exact seven F4 terminal markers without a service warning or alert. The service time and memory fields remain opaque metadata; they are not process-scoped user/system CPU or peak RSS. The basis polynomials and a Magma-linked rational point witness were not obtained, so the status is `sat_basis_certificate_unverified_model`.

This is the first retained Magma F4 terminal at `n=41` in the campaign. It strengthens solver-coverage evidence at a literature-scale PDP degree, but it is one source-equivalent planted instance rather than five process-metered replicates. It does not execute `n=59`, complete the Magma matrix, establish end-to-end index-calculus scaling, or change the non-SOTA conclusion.
