# Rotated m5/m6 solver admission: blocked at complete-instance export

**Decision: zero solver arms admitted from the pinned interfaces.** The fixed
n13-m5/n19-m6 corpus is usable: each rung has four planted and four
independently certified projected-negative subgroup Q targets, with all four
`Q+T` point branches recorded. The current examined S3-chain, binary
resultant, SAT and WDSat interfaces do not yet export a *complete, matched*
rotated point-PDP instance for these labels. No solver was run on a PDP input;
there is no SAT/UNSAT solver verdict, timing ratio, full-DLP S or rho claim.
This is an **interface-admission blocker for the hash-pinned paths below**, not
a claim that a complete exporter is impossible or absent everywhere in the
repository.

The preregistration is in [PROTOCOL.md](PROTOCOL.md); the initial draft
`c901d44` passed hash-only CI before the first static child. That child's
binary inventory covered PATH only. We retained its complete raw receipt,
then amended and re-froze the source to include the existing #764 WDSat
fixture outside PATH. Revised draft head `7c7bee1` passed both hash-only
workflow jobs and the workflow parser before the final child at
`14:24:26.381929–14:24:26.686720 UTC`. The final `FROZEN.json` SHA-256 is
`e8159d2b21ed2da4b8591d792ed1c3c297a38933addb66a86c912ef2070e374b`;
the source/evidence inputs are SHA-pinned in
[INPUT.json](INPUT.json). A fresh archive-only replay passed. No measured
solver campaign or trial was hidden by this inventory amendment.

| Audited arm/interface | Exact frozen corpus/control | Current admission limit | Solver status |
|:--|:--|:--|:--|
| n13-m5 rotated, 8 Q | 3,125 point tuples; four planted/four certified negatives; 32 `Q+T` branches. Affine S3 paths per Q: `2,2,2,2,0,0,0,0`. | Raw chain layout is 49 bits, within the generic 64-variable mask, but `build_decomposition_system` uses one shared basis for all slots and exports no O/inverse or rational-lift branch. No frozen rotated input or complete model map exists. | **Not admitted** |
| n19-m6 rotated, 8 Q | 117,649 point tuples; four planted/four certified negatives; 32 `Q+T` branches. Affine paths: `3,1,4,2,0,0,0,0`. | Raw chain layout is 88 bits before auxiliaries, above the generic 64-variable mask; the same rotated-domain and exceptional-branch export is absent. | **Not admitted** |
| Direct binary resultant | `binary_semaev.rs` exposes S3 and S4, with no frozen S6/S7 Boolean or solver file for these 16 targets. | An algebraic zero needs rational two-sign factor lift and full Q+T replay. | **Not admitted** |
| SAT/XOR, WDSat, F4/FES, msolve | Solver adapters and runnable binaries exist for other systems; the S5 SAT example has four factors. | A solver cannot repair a missing complete input. FES also requires a verified quadraticization of S3 or higher-degree rows; msolve lacks this corpus's pinned GF(2) exporter/parser. | **No comparable arm** |

The #777 finite-fibre theorem establishes that an affine S3 root on complete
rational, sign-paired factor fibres corresponds to a rational signed sum, and
that affine chains admit consistent signs. This removes a possible
*affine-root soundness* objection, but does not encode prefixes that equal O.
The #774 n13 exact branch `Q+O=(7256,3272)` has three true x-mask witnesses;
affine S3 finds two and misses exceptional-only mask `[0,0,0,2,1]` whose
first two factor points are both `(0,1)`. That Q remains affine-SAT on this
specific corpus. The missing mask proves model-set incompleteness and explains
why an affine-only UNSAT rule cannot be promoted to general projected PDP
UNSAT. The four negative Q per rung are certified by the *independent point
census*, not by a SAT or Gröbner refutation. Every one of their four torsion
cosets is empty in the inherited archive and matched the #774 replay.

Local binary inventory was charged in the final static child: CryptoMiniSat
5.14.7, msolve 0.9.5, Kissat 4.0.4 and CaDiCaL 3.0.1 were on PATH; MiniSat
and `wdsat_solver` were not. The #764 WDSat fixture did exist outside PATH at
`/private/tmp/kic-wdsat-fixture-20260925/wdsat_solver`, SHA-256
`fe874e42bccda6125b6588cb2840a58c3a31d3ad13b48c1a2047f7cc098ec249`,
matching its earlier committed receipt. All on-PATH binary hashes, paths,
version output and exit statuses are in `evidence/final/result.json`; a local
binary's availability does not establish an admitted problem encoding.

The final static child exited zero in **0.186 s wall / 0.045 s CPU**, with
25,886,720-byte self RSS and a 25,935,872-byte child high-water bound, below
the frozen 30-second/256-MiB caps. Its raw `result.json` SHA-256 is
`7d1acc98e428af3674e3f3d875f223a825dbdd66e0b8a8702e0bb93c6857dc8b`;
receipt SHA-256 is
`7e2eeb25931d8bb4aca2f2289e3fd4de1445c91b8d945f4ee69dbff3a646ce17`.
The archive replay recomputed the 16-Q/64-branch counts and bounded-interface
facts without rerunning any heavy oracle. [Raw evidence and replay](evidence/README.md)
retain both the original and corrected attempts. Classification: **accounting /
solver-admission blocker**; S and matched-rho ratios remain unset.

The next experiment must implement one complete rotated encoding on the same
16 Q targets, beginning with an n13 chain that includes explicit infinity,
inverse and O-target branches, both rational factor signs, all four torsion
translates, and independently decoded point witnesses. Then handle n19's
88-bit pre-auxiliary chain through a multiword or direct circuit/export path.
Freeze each engine's CNF/XOR/ANF bytes and model/exit parser, verify its full
projected model set against the point oracle including the #774 exceptional
mask, and only then apply the 15-second/2-GiB cold pilot caps from
[#763](../ROTATED_M5_M6_SOLVER_ADMISSION_20260925.md). Compare positive
first-valid-witness and negative complete-refutation workloads separately,
with failed/censored runs retained. This gate neither prices n131 solver cost
nor changes its factor-base coverage and rho boundary.
