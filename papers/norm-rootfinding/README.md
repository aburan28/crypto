# Exact Norm-Based Search for Constrained Roots over Finite Fields

**Research draft; anonymized; not peer-reviewed; not submitted to a journal,
conference, arXiv, or ePrint; not cleared for submission.**

Change class: engineering and mathematical exposition of existing versioned
components. No priority claim or measured F4/F5/F6 speedup. This is separate
from `../rootfinding-decomposition/` and does not import its measurements.

## Claims and boundaries

The manuscript contains readable proofs of the repeated-root-safe support
predicate, the characteristic-two norm-chart parametrization, unique extension
under bounded affine contraction, and source-permutation restoration.
Its implemented enumeration contract is **norm-primary models**, not original
S3-chain roots. Zero/repeated sources remain. The z=0 target branch, chain
auxiliary reconstruction and ghost/exception closure remain open.

The recorded GF(2^83) profile reduces native internal coordinates
58,176 -> 45,427 (21.9%). The following bounded hard search returns UNKNOWN,
zero verified witnesses, null exact norm count and null original root count.
These are diagnostics, not runtime acceleration. Odd-field tests cover
supplied quintics, including the five NIST prime fields, rather than a completed
prime-field multivariate solver. The full qualified speed gate remains unmet.

## Build and reproduce

Run `make check` here; it uses `../paper.mk` and requires pdflatex, bibtex,
pdfinfo and standard LaTeX packages. Rebuild presentation figures with
`python figures/render_figures.py` (matplotlib; presentation only).
The compiled PDF is `norm_rootfinding.pdf`.

`data/evidence_map.json` maps claims to immutable source paths, byte counts,
and SHA-256 digests at commit 4619bb97f63d960d40ffce127154b4f54195135d.
Four raw records are copied unchanged. `data/fixture_metadata.json` supplies
the binary norm parameters; the full bulk original-chain input was excluded
from publication and remains a concrete prerequisite for its byte-identical
replay. Frozen native source/input manifests remain authoritative.
No solver experiments or benchmarks were rerun to generate this paper.

## Verification status and missing work before submission

- Readable central coefficient proof is included; native literal-norm controls
  are cited. No theorem here is Lean-certified. `lean` and `conductor` were
  unavailable locally. See `PROOF_MAP.md` for exact outstanding declarations.
- The current identity.certificate/v1 evaluator fixes an odd prime and cannot
  certify characteristic-two/conditional claims faithfully. No false
  certificate is issued. A characteristic-aware certificate adapter remains.
- Complete the novelty audit, especially WDSat (Trimoska/Ionica/Dequen),
  divisor/function-field parametrizations, and later algebraic SAT work.
  Bibliographic verification and abstract review do not equal a full-paper audit.
- Obtain independent mathematical review and independent evidence replay.
  This draft's author performed self-review, not independent peer review.
- Prove coverage of the original system, all auxiliary fibers and exceptional
  branches; implement an independently checked exhaustion certificate.
- Prove and implement an odd-characteristic search chart and nonbinary domain
  contracts. Supplied-polynomial roots do not fill this obligation.
- Complete source-identified genuine F4/F5/F6 and strongest relevant external
  comparisons on identical full all-roots outputs, including field coverage,
  20 hard instances, >=3 cold runs, qualified isolation, A/A <=10%,
  >=2x geometric floor and <=1.10 worst slowdown. Requested higher speed
  targets remain open; no stage count or timeout may qualify release.
- Next experiment: sorted-source comparator with complete unique permutation
  restoration; then whole-cell rejection and conditional factorization.
  Freeze changed sources and budgets; do not rerun the identical failed search.

The earlier omitted-term transfer formula is corrected explicitly in the paper,
with the GF(4) counterexample and preserved erratum linked. No cryptographic
key-recovery, security break, or asymptotic ECDLP claim is made.
