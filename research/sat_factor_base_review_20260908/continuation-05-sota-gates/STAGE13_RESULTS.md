# Stage 13: repeated target-matched PDP scaling panel

The frozen panel ran five unchanged seeds for four planted point-decomposition cells: `n=31,ell=5,m=3` standard, the selected nondegenerate `n=31` GGMP kernel, `n=41,ell=5,m=3` standard, and `n=59,ell=9,m=3` standard. All 20 tasks ran sequentially from clean commit `67cf7f5` under a whole-cell process meter. Every solver arm was bound to one target, public factor-base predicate, regenerated source system, producer BLAKE3 identity, and pre/post SHA-256 export custody.

Native SAT and direct MITM ran in isolated fresh processes. Native SAT blocked nonlifting algebraic x-tuples and accepted SAT only after a rational factor-point sum matched the target. WDSat and CryptoMiniSat SAT models passed separate metered point-witness validation. Magma inputs were generated with direct sparse F4, one thread and GPU disabled, but no Magma executable was available.

The table reports medians across five source-bound replicates. Solver core time excludes the separately listed WDSat/CryptoMiniSat witness-validation process; the whole-cell totals include every sequential descendant and Python orchestration.

| Cell | Backend | Terminal outcomes | Median conflicts | Median core-s | Median wall-s | Median peak MiB | Median witness-validator core-s |
|:--|:--|:--|--:|--:|--:|--:|--:|
| n31 standard | native XOR | 5 SAT | 1,855 | 0.061549 | 0.066491 | 5.64 | - |
| n31 standard | WDSat | 5 SAT | 2,212 | 0.015315 | 0.399260 | 3.48 | 0.012445 |
| n31 standard | CryptoMiniSat | 5 SAT | 15,179 | 0.194419 | 0.245184 | 7.80 | 0.009761 |
| n31 standard | direct MITM | 5 SAT | - | 0.013175 | 0.017678 | 4.09 | - |
| n31 GGMP | native XOR | 5 unknown at 100k | 100,000 | 2.460058 | 2.913967 | 23.53 | - |
| n31 GGMP | WDSat | 5 timeouts | - | 98.595602 | 120.007911 | 7.28 | - |
| n31 GGMP | CryptoMiniSat | 5 SAT | 620,648 | 13.141137 | 13.904087 | 54.97 | 0.007714 |
| n31 GGMP | direct MITM | 5 SAT | - | 0.022424 | 0.030054 | 3.56 | - |
| n41 standard | native XOR | 5 SAT | 6,252 | 0.203938 | 0.253097 | 7.52 | - |
| n41 standard | WDSat | 5 SAT | 2,520 | 0.027199 | 0.559454 | 3.67 | 0.012174 |
| n41 standard | CryptoMiniSat | 5 SAT | 16,992 | 0.231797 | 0.257680 | 8.98 | 0.013035 |
| n41 standard | direct MITM | 5 SAT | - | 0.020080 | 0.025243 | 4.47 | - |
| n59 standard | native XOR | 5 unknown at 100k | 100,000 | 11.989420 | 14.357174 | 157.14 | - |
| n59 standard | WDSat | 5 timeouts | - | 98.396824 | 120.016037 | 13.58 | - |
| n59 standard | CryptoMiniSat | 5 timeouts | - | 96.517291 | 120.013840 | 129.17 | - |
| n59 standard | direct MITM | 5 SAT | - | 3.857910 | 4.371598 | 44.24 | - |

One original `n=59` WDSat process aborted before search because its generated static buffer was 32,000 entries for 32,263 global XOR atoms. The immutable failure remains in the panel. An additive same-seed correction sized the buffer to 32,264, reproduced the exact target and source/export identities, and reached the 120-second watchdog. The corrected distribution therefore contains five WDSat timeouts. The correction process consumed 114.916149 core-seconds, 148.127927 wall-seconds, and 160.88 MiB peak RSS; its repeated native/MITM/build costs are retained.

The original 20 task envelopes consumed 1,466.977905 core-seconds and 1,946.021003 sequential wall-seconds. Charging the one-time twelve-candidate GGMP discovery and the additive WDSat correction brings the recorded campaign total to 1,582.114824 core-seconds and 2,094.978436 sequential wall-seconds, with 242.95 MiB maximum RSS.

The panel gives repeated scaling evidence for planted PDP instances through `n=59`. It does not show a SAT advantage: direct MITM has the lowest median core time in every cell, and all SAT solvers are inconclusive at `n=59` under their caps. It is not a balanced random SAT/UNSAT panel, an end-to-end index-calculus run, or a SOTA result. The full solver-matrix gate remains open because Magma was unavailable.
