# Frozen protocol: exact n19 rotated-base portfolio support

Status: preregistered before the first five-base or higher-order union outcome.
This is an archive analysis of the independently replayed full subgroup census
in [PR #767](https://github.com/aburan28/crypto/pull/767) and four equal-size
normal-generator arms in [PR #769](https://github.com/aburan28/crypto/pull/769).
Those prior reports already disclosed single-base sizes and beta-3 pairwise
overlaps. This protocol fixes the further calculations, code, and inputs before
inspecting the five-base result. No base is selected or replaced after outcome.

## Question, population, and measures

Use the exact n19 Koblitz group with prime subgroup order `q=130873`,
`H=(385982,301867)`, polynomial `x^19+x^5+x^2+x+1`, m=6, d=2,
and the ordered normal generators `[3,338435,303097,464276,42605]`.
Each factor-base arm has seven physical points per slot and three nonzero
signed cofactor-projected F0 columns plus O. The target population is **all**
`q` projected subgroup points, not the eight fixed labelled Q targets alone.
Membership means a positive count in that arm's archived complete projected
six-sum histogram; the same point can have different multiplicities by base.

Compute the exact 32 membership-pattern frequencies (including misses for all
five), all 31 nonempty subset union sizes, and in the fixed order above each
prefix's newly supported points and remaining misses. Also compute the mean
number of sequential membership-oracle probes in that order, assuming a free,
perfect oracle that stops on first supported base and tests all five on misses.
Its exact numerator is `sum_{j=0}^4 misses(prefix_j)`, with the empty prefix
having `q` misses. This is a **query-count diagnostic** only: it omits
membership-test runtime, model search, relation collection, linear algebra,
new factor-base setup, and logarithm verification.

The support-only follow-up gate is five-base union at least `ceil(0.90*q)`.
If it passes, choose the smallest subset crossing this threshold, tie-breaking
by the fixed beta order, solely as a candidate for a later charged experiment.
If it fails, do not prioritize this five-base portfolio. Passing cannot imply
an attack speedup: each new base may require its own independent relation rank
and logarithms, and all costs must be charged before comparison with rho.
The selected subset must be tested on new held-out problems or a larger rung
before any transfer statement about n131.

## Inputs, independent replay, and bounds

The primary analysis reads the five archived projected point histograms and
recomputes all set operations without group arithmetic. It verifies each
point occurs once in a histogram, its positive multiplicity, the per-arm
`7^6` tuple total, and known single/pair support sizes. The independent
verifier uses the separately implemented bit-serial/Fermat field and group
law from #767 to enumerate `k*[4]H`, `0<=k<q`. For the four #769 arms, it
reads the independent `target_counts.u32le` arrays and checks each k count
against the point histogram before recomputing all membership patterns and
subset unions. The beta-3 arm has no integer array, so it maps that arm's
archived point histogram through the same verified group enumeration. It
compares the complete outcome JSON, not a sample. Both inputs are pinned by
SHA-256 and both scripts are frozen before the first portfolio output.

Run analysis and verifier as separate cold children with 120-second wall and
512-MiB RSS caps each. Retain UTC intervals, exit status, stdout/stderr,
source/input/result hashes, and any failure or censor receipts. An over-cap
or failed verifier leaves the portfolio decision unknown. CI first checks only
the frozen hashes and source syntax; after measurement it checks archived
results and reruns the independent verifier. No solver, relation rank, ECDLP,
or matched-rho speed ratio is measured in this experiment.
