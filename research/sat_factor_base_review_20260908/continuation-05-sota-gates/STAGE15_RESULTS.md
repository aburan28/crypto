# Stage 15: bounded official Magma Calculator evidence

Two of the twenty frozen Stage 13 Magma inputs fit the official calculator's observed 50,000-byte limit: the `n=31,ell=5,m=3` standard instances for seeds `2026091301` and `2026091303`. The other eighteen frozen inputs are 61,062 to 447,248 bytes, fall outside that cap, and have no retained calculator responses. The retained XML identifies Magma V2.29-10 and an enforced 60-second limit.

The byte-exact archived scripts both continued through F4 and emitted complete terminal markers, but the calculator rejected their leading global `SetNthreads(1)` call as an illegal sandbox operation. The strict Stage 13 parser therefore does not admit those exact responses as clean solver terminals. Additive calculator-specific inputs removed exactly that one line, retained `SetGPU(false)` and `Nthreads := 1` on the direct sparse `GroebnerBasis` call, and left the Boolean ring and equation sequence byte-identical.

| Seed | Frozen input bytes / SHA-256 | Adapted input bytes / SHA-256 | F4 result | Step degrees | Basis size | Internal CPU-s | Internal wall-s | Service time-s | Service memory field |
|:--|:--|:--|:--|:--|--:|--:|--:|--:|:--|
| 2026091301 | 33,822 / `fd82e5b47caf1904bd9e7cc5e6c8ae126644014868ea23a396b8f2bb5f3bacfd` | 33,806 / `1de047cafe10b89eaeef254b291dcc61b991a6edafc1c64d1d72bdd08a115101` | proper ideal / SAT certificate | `[2,3,3,3,3]` | 102 | 0.440 | 0.450 | 0.480 | `32.09MB` |
| 2026091303 | 26,834 / `93dbea20524f854176e4bfc3ef68080b3a8e5b9df63f206016a79d8c7e2a9689` | 26,818 / `95103be43ce34379aa2a88afd360a74dc72e98b30b50c35286b732d24b6723d4` | proper ideal / SAT certificate | `[2,3,3,3,3]` | 102 | 0.420 | 0.420 | 0.450 | `32.09MB` |

The service's time and memory headers are retained as opaque service-reported fields. They are not process-scoped user/system CPU or peak RSS and cannot close the resource-accounting gate. The internal CPU and wall markers surround only `GroebnerBasis`.

Separate witness attempts did not produce assignments. Magma's `SAT` wrapper failed because the calculator sandbox disallowed `GetTempDir`; a documented variety attempt and an attempted full-basis export both reached the 75-second client watchdog without a response. Those are operational failures and assert nothing about point decomposition, variety size, or basis-export cost. Consequently both clean F4 results remain `sat_basis_certificate_unverified_model`, not validated SAT point decompositions.

The fail-closed verifier checks the frozen Stage 13 SHA-256 custody, the exact one-line calculator adaptation, raw XML hashes and headers, agreement between exact and adapted F4 terminal identity fields (status, step degrees, basis size, algorithm and controls), the failed witness boundary, and the timeout receipts. It regenerates `stage-15-summary.json` byte-for-byte.

This is useful partial Magma evidence: two target-matched systems now have calculator-reported direct sparse F4 proper-ideal terminals from an official external service. Their reported status, step degrees and basis size agree; the basis polynomials themselves were not retained because the export attempt timed out. This is not the twenty-instance Magma matrix, lacks required process metrics and rational point witnesses, is not unaffiliated methodology or novelty review, and does not change the non-SOTA conclusion.

References: the official [Magma Calculator](https://magma.maths.usyd.edu.au/calc/) states its input/runtime limits and version; the official [Gröbner-basis documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/1314) defines direct Faugère F4, `Dense`, `Nthreads`, Boolean polynomial rings, and the returned F4 step-degree sequence; the official [SAT interface documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/1318) describes source-model extraction separately from Gröbner-basis construction.
