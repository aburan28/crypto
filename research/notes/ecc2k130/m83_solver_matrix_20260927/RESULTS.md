# m=83 Koblitz index-calculus solver evaluation, 2026-09-27

The frozen runs use `E0: y²+xy=x³+1` over `GF(2^83)` modulo
`z^83+z^45+z²+z+1`, with subgroup order
`r=2417851639230796216685689`, cofactor 4 and Frobenius order 83.
This evaluates decomposition oracles and a tiny-factor-base relation search.
**No discrete logarithm was recovered.** End-to-end costs, rho comparison,
speedup, and transfer to ECC2K-130 remain unestablished.

## Matched fixed-phase quotient workload

Each engine received the same six-variable Semaev-S4 Boolean system for each
of 12 inputs: two seeds, two planted targets and four independently sampled
natural targets per seed. The factor-base payload has two bits; phases are
fixed to `(0,1,2)`. Median wall seconds measure *solver setup, solving and
complete algebraic model extraction only*. They exclude equation build,
lifting, factor-base setup and any DLP work. FES enumerates 64 assignments.
The Python F4-style engine uses batched squarefree Macaulay elimination;
the Python F5B engine uses signatures plus a completion check. These are
explicitly prototypes, not timings of Magma F4/F5 or msolve.

| Engine | Completed / planned | Median solver-stage wall s | Natural relations / 8 | Correct versus FES |
|:--|--:|--:|--:|:--|
| FES, exhaustive reference | 12/12 | 0.000335 | 0 | independent reference |
| CryptoMiniSat 5.16, CNF + native XOR | 12/12 | 0.0281 | 0 | 12/12 |
| Boolean F5B prototype | 12/12 | 0.1497 | 0 | 12/12 |
| Boolean Macaulay/F4-style prototype | 12/12 | 0.1710 | 0 | 12/12 |
| Z3 4.15.3, Boolean And/Xor | 12/12 | 0.2121 | 0 | 12/12 |

CryptoMiniSat had the shortest measured SAT-stage wall time in this restricted
workload. This does not rank complete m=83 index-calculus algorithms. Two
algebraic assignments across these 12 distinct inputs failed group lifting;
all engines that completed returned the same full root sets and rejected the
same false positives. The planted positive controls verified. Natural targets
produced no group relations. No operation-unit conversion between SAT,
matrix rows, F5 pairs and FES assignments was measured, so their cross-engine
*operation* ratios are null; wall times are diagnostics.

## Separate normal-subspace payload sweep

The archived cyclotomic factor-base constructor requires `s | (83-1)`.
Attempts at `s=3,4` failed at setup and are preserved in the preflight
receipt. We therefore preregistered a separate nested normal-basis subspace
factor base for `s=2,3,4`. Its `s=2` inputs are *not* the quotient inputs
above. At each size and seed the grid contains one planted control and one
natural target per engine. A setup failure affects the same target for all
engines. FES tests 64, 512 and 4096 assignments, respectively.

| Payload bits | Engine | Complete / 4 | Timeout / 4 | Setup failed / 4 | Median completed solver-stage wall s |
|--:|:--|--:|--:|--:|--:|
| 2 | FES | 2 | 0 | 2 | 0.000311 |
| 2 | CryptoMiniSat XOR | 2 | 0 | 2 | 0.0319 |
| 2 | Z3 | 2 | 0 | 2 | 0.2155 |
| 2 | F5B | 2 | 0 | 2 | 0.1215 |
| 2 | F4-style | 2 | 0 | 2 | 0.1957 |
| 3 | FES | 2 | 0 | 2 | 0.0100 |
| 3 | CryptoMiniSat XOR | 2 | 0 | 2 | 0.0311 |
| 3 | Z3 | 2 | 0 | 2 | 0.7166 |
| 3 | F5B | 0 | 2 | 2 | — |
| 3 | F4-style | 0 | 2 | 2 | — |
| 4 | FES | 3 | 0 | 1 | 0.3525 |
| 4 | CryptoMiniSat XOR | 3 | 0 | 1 | 0.2147 |
| 4 | Z3 | 0 | 3 | 1 | — |
| 4 | F5B | 0 | 3 | 1 | — |
| 4 | F4-style | 0 | 3 | 1 | — |

Every completed solver result matched the exhaustive algebraic root set and
independent group checks. No natural target yielded a verified relation.
Most planted controls could not be generated from the frozen plain subspaces:
the first seed had no usable subgroup lifts through `s=4`; the second seed
had a planted positive only at `s=4`. This limits inference from the sweep;
we preserved those cells and did not replace the seeds. CryptoMiniSat's
binding exposes no seed setting; it ran one thread with its version's default.

## Phase-choice encodings and complete orbit control

The archived phase-choice probe allows three phases per point. It uses two S3
equations with a field-valued intermediate, so it is a **different input
representation** from the six-variable direct S4 comparisons. On both seeds,
the inline arm had 95 variables, roughly 192–196 thousand Boolean terms,
about 7 seconds of SAT setup and a five-second Z3 `unknown/timeout`.
The implicit arm had 347 variables, about 3.49 million Boolean terms and
roughly 22 seconds of equation building. Both implicit runs hit the
25-second SAT-setup budget, after peak memory about 920 MB. FES, F4-style
and F5B require an exponential squarefree layout at 95 or 347 variables:
their preflight is `unsupported_memory`. WDSat, msolve and Sage executables
were unavailable; no numbers are imputed to them.

For the original quotient factor bases, exhaustive meet-in-the-middle over
**all 83 Frobenius phases and both signs** gave:

| Seed | Signed orbit | Pair sums tried | Group additions, all charged stages | Planted controls | Natural relations / 4 | Natural rank |
|--:|--:|--:|--:|--:|--:|--:|
| 260938 | 332 | 55,278 | 94,592 | 2 verified | 0 | 0 |
| 260939 | 498 | 124,251 | 181,011 | 2 verified | 0 | 0 |

The counts include orbit verification, pair construction, target search and
verification, and setup group additions. They do not convert omitted field
setup operations into group additions. Each planted target had six verified
triples under sign/permutation symmetry; one planted target per seed had six
tautological presentations. All eight natural targets had zero full-orbit
relations. There is no relation matrix of useful rank and no final scalar
recovery. A measured rho operation reference and `S=total/sqrt(r)` therefore
remain null.

## Reproduction and decision

`results/validated_summary.json` is checked by `summarize.py` against
every available FES receipt, including hashes, full algebraic root sets and
group verification results. The archive named in `results/README.md`
contains every subprocess wrapper, source hash manifest, preflight failure,
timeout, and full-orbit transcript. `run_v1.py` exactly matches the
original run's recorded source hash; the current `run.py` implements the
subsequently preregistered extension. Resource limits were 1 GiB virtual
address space, 25 seconds equation construction/setup where applicable,
five seconds solving and 70 seconds per child.

The bounded result is: native-XOR SAT is the fastest measured SAT *stage* for
these fixed-phase diagnostic inputs; tiny FES is faster at six variables;
Python F4/F5 prototypes lose completion at the larger payloads. Neither
those rankings nor the zero-relation tiny-base sample establishes an
end-to-end m=83 IC improvement. The required high-fidelity ECC2K-130 gate
remains open pending a rank-complete, verified-DLP, full-cost comparison on
matched m=83 workloads and a measured rho reference.
