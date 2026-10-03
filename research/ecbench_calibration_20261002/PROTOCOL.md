# ecbench calibration: does the harness charge each method what theory says?

**Preregistered 2026-10-02, before any calibration session ran.** This
file and the two specs in [`specs/`](specs) were committed before the
first session. Results go in [`README.md`](README.md) and
[`sessions/`](sessions), never back into this file.

## Question

`ecbench` charges every generic ECDLP method through one counted ledger
(`docs/ecbench/README.md` §7). If that ledger is right, each method's
mean `S = total group operations / √r` must approach the constant its
analysis gives, from above, as `r` grows. If it does not, the harness
(or a method's implementation) is wrong, and no comparison made with it
can be trusted.

This is a **validation of the measurement instrument**, not an attack
claim. It changes no method, and no ratio here is a speedup.

## Boundaries (derived, written before measuring)

The floor of each curve is `√(π/2A)` in `S`, with `A = 2` on the prime
curves and `A = 2m` on the Koblitz curves. The expectation of each method
(`methods::expected_s`) is:

| method | expected `S` (leading order) | basis |
|---|---|---|
| `rho.plain`, `rho.frozen_reference` | `√(π/2) ≈ 1.2533` | birthday bound on `r` points |
| `rho.negation` | `√(π/4) ≈ 0.8862` | classes `{P, −P}` |
| `rho.signed_frobenius` | `√(π/4m)` | classes `{±φ^t P}`, `A = 2m` |
| `bsgs.textbook` | `1.5` | `m = ⌈√r⌉` baby steps plus a mean of `m/2` giant steps |
| `bsgs.interleaved` | `4/3 ≈ 1.3333` | `E[2·max(i, j)]` for uniform `i, j < m` |
| `bsgs.negation` | `1.0` | `√r/2` baby steps plus a mean of `√r/2` giant steps of stride `2m+1` |
| `kangaroo.vow` | `2.0` | van Oorschot–Wiener, one tame and one wild, mean jump `√r/2` |

Each omits `O(log r)` set-up (jump tables, the giant stride's scalar
multiplication), so the measured `S` should exceed it at small `r` and
approach it as `r` grows.

## Frozen inputs

- Specs: [`specs/prime.json`](specs/prime.json) (six prime-order curves,
  `find_prime_order_curve(N, 59297)` for `N = 16, 18, …, 26` as
  `docs/curves/registry.json` records them, given explicitly) and
  [`specs/koblitz.json`](specs/koblitz.json) (six Koblitz curves,
  `r = 2^16 … 2^39`). Every slug is registered, with an EC1 alias.
- 8 independent single-target workloads per curve (`target_seed`
  20261002), 3 rounds each (algorithm seed 20261002), no warm-up (counts
  do not depend on warm-up), arms interleaved by `alternate`.
- The `ecbench` binary built from the commit that adds this file, release
  profile. Its SHA-256 is recorded in every session.

## Metric and statistics

For each (method, curve): the mean `S` over its 24 verified runs, with a
95 % two-stage bootstrap interval (8 workloads resampled, then runs within
each; 4 000 resamples, seed 20261002), as `ecbench table` prints it. The
quantity tested is `S / theory`.

## Predictions and falsification

At the **largest curve of each family** (prime `r ≈ 2^25.4`, Koblitz
`r = 2^39`):

1. **BSGS (exact analysis).** For `bsgs.textbook`, `bsgs.interleaved` and
   `bsgs.negation`, the theory constant lies inside the 95 % interval of
   the mean `S`.
2. **Rho (heuristic random-walk analysis).** For `rho.plain`,
   `rho.negation` and `rho.signed_frobenius`, `S / theory` lies in
   `[0.85, 1.35]`. An r-adding walk is a few percent less random than a
   random map, and the negation and Frobenius walks add fruitless-cycle
   handling, so `S / theory` slightly above 1 is expected.
3. **Kangaroo (heuristic).** `S / theory` lies in `[0.8, 1.5]`.
4. **Correctness.** Every execution verifies (`[k]G = Q` checked by the
   runner), and `ecbench verify --replay 12` reproduces every replayed run
   exactly, on this host and in CI on Linux x86-64.

The approach to the constant across sizes is reported, not tested. With 8
targets the sampling spread of a BSGS mean (about ±0.1 in `S`) exceeds the
`O(log r / √r)` set-up term such a test would resolve.

`rho.frozen_reference` is reported and not tested: it is the historical
"before" walk, kept to show what the tuned walks changed.

**Stop and investigate** if any execution fails to verify, any replay
differs, or prediction 1 fails. A failure of prediction 1 is an accounting
defect in the harness or the method: it is fixed in a separate commit,
classified as *accounting* (AGENTS.md §3), and the calibration reruns as a
new session. The failed session is kept.

**Inadmissible:** changing curves, targets, seeds, the expected constants
or the tolerance bands after seeing results; dropping or rerunning
individual executions; pooling the two families.

## Host and isolation

The sessions run on the author's Apple M4 Pro (macOS, P and E cores, no
affinity control), so every run is **L0**. That is sufficient for
operation counts, which do not depend on the host. Wall time is recorded
and is not a result. The CI job replays the committed sessions on Linux
x86-64, which is the cross-host check of the counts.
