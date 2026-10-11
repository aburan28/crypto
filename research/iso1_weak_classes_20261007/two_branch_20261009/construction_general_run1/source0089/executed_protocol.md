# Two-branch classification and representative construction

Registered 9 October 2026 before executing this round. The approved target is
an exact criterion for both linear and quadratic branches of the genus-3
hyperelliptic family, an explicit route to a representative, and construction
costs from independently sampled classes. This round preserves the existing
192–252-bit field requirement and the larger-prime census as separate tracks.

The family is `y²=h(x)(x-alpha)(x-alpha^q)` over `K=F_(q³)`, `q=p²`,
`alpha` outside `F_q`, with squarefree `h` of degree one or two. Split
quadratics are to be reduced to the linear branch by an `F_q` Möbius map;
nonsplit quadratics are normalized to `h=x²-d` for a fixed base-field
nonsquare `d`, allowing both quadratic twists. These normalizations require
proof and executable checks before an exhaustive family label is assigned.

- Enumerate the cubic torus of order `p⁴+p²+1` modulo inversion and absolute
  Frobenius, and the proposed quadratic torus of order `p⁴-p²+1` modulo
  absolute Frobenius. Check the quadratic cross-ratio formula and field
  projection for every retained parameter.
- Begin at `p=7`, comparing cubic support with the existing exact all-x
  census. Check the quadratic parametrization against all `p⁶-p²` possible
  normalized alpha values, using invariant comparison and separate count
  controls. Expand complete two-torus enumerations to feasible larger primes.
- Compare the combined class labels against independently sampled source
  curves. The first validation workload reuses the 64 previously frozen
  independent `p=7` sources, not the forced-positive fixtures. A source
  population is model-uniform, not trace-uniform; repeated classes stay rows.
- A constructive route must consist of explicit rational isogeny maps and a
  verified final model conversion. Count every attempted edge and visited
  vertex. Record setup, construction, and verification separately. Finite
  caps produce unresolved rows, never certified class exclusions.
- Evaluate a direct quadratic membership test and broaden the existing
  192–252-bit search to quadratic endpoints if correctness validation passes.
  Keep all frozen prior evidence intact and retain every new capped row.
- Execute arithmetic in compiled PARI/GP, with native Rust orchestration,
  checking, and identity records. Preserve script, binary and protocol hashes,
  raw stdout/stderr, exact field moduli, inputs, and all process statuses.
  Resource wall times on this shared host are construction diagnostics;
  controlled timing comparisons require the host-isolation gate.
- Publish the theorem's exact domain, complete proof, identity certificates,
  dated report, explanatory diagram, quantitative graph and reviewed PDF in
  the existing PR. General large-field classification stays partial whenever
  the exact decision procedure exceeds resources or remains unimplemented.

The cubic-only necessary condition is not an exclusion criterion for the
combined family. The previous high-depth cubic zeros retain their original
meaning and may have quadratic representatives.
