#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.5pt)
#set par(justify: true)
#show table: set text(size: 8.5pt)

= lit2026 W12: flagged ECDLP claims

*Zero-minor guess-and-determine (arXiv 2607.09814) against rho.negation; Hu ePrint 2025/015.*

2026-10-11 · audit workstream W12 · run `run-20261012T0009Z` · solver commit `c70d57ec7` · L0 host · counted field operations.

== Verdict

#table(columns: (auto, 1fr),
  [*Zero-minor method* (Mahalanobis; 1703.07544 → 2310.04132 → 2607.09814)],
  [*Negative.* Implemented as specified; it solved 48/48 planted targets on the registered 16–26-bit prime ladder.

  - *Scaling* in counted $F_p$ ops, against $log_2 r$:
    - search phases: slope 0.905 (95% bootstrap 0.769–1.021);
    - all phases: 0.783 (0.669–0.879);
    - candidate minors examined: 0.929 (0.793–1.048), with a median of 0.37–2.1 $times r$ per curve.
  - *Constant:* 260× (215–319) to 6715× (2263–15336) the cost of `rho.negation` per target.
  - *Falsifier met.* The slope is ≥ 0.5 and the constant is far above rho.
  - The paper's own Table 2 says it is "no better than exhaustive search".],
  [*Hu ePrint 2025/015*],
  [*Blocked.*
  - The full text is behind an IACR Cloudflare challenge, which was not bypassed.
  - The abstract leaves "IICA" undefined, as well as "integer up/down".
  - The abstract concerns the DLP in $F_p^*$, not the ECDLP.
  - Nothing could be implemented. The original requirement stays unresolved.
  - The companion ePrint 2024/1923 (abstract) says no subexponential prime-field ECDLP method is known.],
)

== The method, in one paragraph

- *Kernel.* Pick $ell = 3n'$ points $n_i P$ and $ell$ points $-n_j Q$. The matrix $M$ of all degree-$n'$ monomials at these $2 ell$ points has a left kernel $K$ of dimension $ell$ (extra-points theorem).
- *Zero minor ⇒ logarithm.* A maximal zero minor of $K$ is a kernel vector with $ell$ zeros. Its support is $3n'$ points on one degree-$n'$ curve, so they sum to $O$, which gives $sum n_i - m sum n_j equiv 0$.
- *Guess.* 2607.09814 guesses $ell' = ell slash 2$ dense columns, transposes and re-reduces.
- *Determine.* For every set $b$ of $ell' - d$ dense columns, take the $d$-dimensional left kernel $T$. Build the $(ell' + d) times d$ signature matrix. Hash the 1-dimensional kernels of all its $(d-1)$-row submatrices; a repeat is a zero $d times d$ minor.
- *Defect.* $d$ is the least value giving success ≥ 0.9 by the paper's Eq. (4): 3 at 16–24 bits, 4 at 26 bits.

== Results

#figure(image("reports/run-20261012T0009Z/scaling.svg", width: 100%), caption: [Run-20261012T0009Z, registered ladder `find_prime_order_curve(b, 59297)`, b = 16…26, 8 targets per curve. (a) Counted field ops per solved target, with OLS fits (95% bootstrap). (b) Candidate $d$-subsets covered, with the line covered = $r$ for reference. (c) Paired per-target ratio of zero-minor to rho.negation (geometric mean, 95% bootstrap).])

#table(columns: 6,
  [bits], [curve], [$r$], [rho.negation median ops], [zero-minor median ops], [ratio zm/rho (95% bootstrap)],
  [16], [`icv1-fp16-t295-8d3c3165`], [53117], [1.92e3], [5.53e5], [260 (215–319)],
  [18], [`icv1-fp18-tm337-d28d5e09`], [212437], [2.98e3], [1.48e6], [455 (326–626)],
  [20], [`icv1-fp20-t727-cd198a38`], [562427], [4.62e3], [2.54e6], [563 (350–883)],
  [22], [`icv1-fp22-tm1385-475dcb5f`], [2435809], [7.93e3], [9.03e6], [785 (349–1633)],
  [24], [`icv1-fp24-t1059-6df599df`], [16482251], [1.67e4], [2.16e7], [1109 (706–1743)],
  [26], [`icv1-fp26-tm1775-7e8fb6df`], [42653927], [2.49e4], [2.71e8], [6715 (2263–15336)],
)

*Identity checks.*
- Every one of the 48 zero-minor and 192 rho runs recovered the planted logarithm. Each was checked by $[m]G = Q$ outside the solver.
- A full replay reproduced 240/240 rows with identical counts.
- The defect sweep at 20 bits ($d = 2, 3, 4$, 8 targets each) gave identical counts under two binaries. Medians: 1.21e7, 2.54e6 and 3.32e6 ops.

*Interpretation (labelled).*
- Each candidate maximal minor is an $ell$-subset of points, which sums to $O$ with probability about $1 slash r$.
- The signature machinery tests candidates in bulk, but it does not raise that probability. So the work is $Theta(r)$ candidates.
- The covered/$r$ medians and the search slope near 1 match this.

== Ambiguities and what was not implemented
+ *Signature.* The paper's RREF signature equals $T c$ up to a scalar; the code uses $T c$.
+ *Representative of a line.* The code normalises the last non-zero entry to 1; the paper normalises the first. Both are canonical.
+ *Collisions with more than $d$ rows* are skipped, as in Algorithm 2: 0 occurred.
+ *Rank-deficient submatrices* are treated as zero minors: 0 occurred.
+ *After a failed guess $a$*, the code redraws $a$. Every run succeeded on its first $a$.
+ *$n'$* is the largest even integer $<= log_2 r$.
+ *The "Pollard rho on penultimate intersections" (§5)* defines no walk, so it was not implemented.
+ *Gray-code and rank-one updates* were not implemented. At 26 bits the $T$ computation is only 7.7% of the ops; the $(d-1)$-kernels are 91%.

== Limits
- *Sizes.* 16–26 bits ($r <= 2^(25.35)$), the largest registered ladder curve. The workstream asked for 20–30 bits; 28–30 needs new registered curves (partial).
- *Seeds.* One seed per zero-minor target.
- *Rho's slope.* rho's fitted field-op slope (0.377) is below 0.5 because of the walk's fixed set-up at these sizes.
- *Wall time.* L0 only, not used for any conclusion.
- *Claims.* No deployed-curve or asymptotic claim is made.

Sources: arXiv 1703.07544v3, 2310.04132v1, 2607.09814v2 (full text); ePrint 2025/015 and 2024/1923 (abstract pages only). Code: `examples/lit2026_w12_zero_minor.rs`. Data: `reports/run-20261012T0009Z/`. Details: `README.md`, `HU_EPRINT_2025_015.md`.
