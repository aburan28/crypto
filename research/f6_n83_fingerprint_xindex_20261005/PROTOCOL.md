# n83 F6 fingerprinted x-index gate

Registered before code changes or candidate measurements. This follows the
rejected flat index in #1428, which cut memory but regressed the dimension-10
query. It is a separate candidate on the retained #1421 parent.

## Hypothesis and exactness

Store a 32-bit array index and a 32-bit fingerprint in each 64-bit
open-address slot. The 4,108,723 full x coordinates remain owned by the
signed-sum array. A probe compares the fingerprint first; if it matches,
it compares the **full** x coordinate in that array before declaring a
hit or duplicate. Therefore fingerprint collisions cannot create false
witnesses or suppress a distinct sum. Preserve insertion order, first
representative per x, the separate infinity slot, sign handling, query
order, portable fallback and independent full-group replay. Grow before
the table reaches 70% occupancy and reject indices beyond `u32`.

The baseline is #1421 head `558cfc1dc23d91c4a38ed5d5a3051d78e045e1da`,
source SHA-256
`feb569917be8b526eee743e988bc01850c3bba32a312f272bd50e42313c3df25`.
Its frozen full/small/planted release binaries are
`75009df1bfe412f00390b1d7e707ba5d4cb64d196a143f7888e6bfe235a15d95`,
`8d2bed511af624b35e5bd48d16b9b42c1920ddaea6ef6163e0c01867a2b6a0e4`,
`4dda20b3331d799036c48bc8362f78ff773c3b78aadfb48c090e48a4ab7987b6`.
Use the same Rust 1.93.1 release flags, physical Apple M4 Pro and
`RAYON_NUM_THREADS=1`. Preserve source, binary and output SHA-256 hashes.

Freeze the registered K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, its exact subgroup, standard
cofactor-projected dimension-8/10/12 bases (258/1,048/4,054 usable
points), and the public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`. The full base has
2,027 signed columns, 8,219,485 unordered pairs and 4,108,723 signed
representatives. The uniform-target four-summand coverage ceiling remains
`4.662e-12`.

## Ordered gates and accounting

1. Run nine focused exact geometry tests and the full-base planted
   `[0,2,4,6]` control through portable and PMULL paths. Every returned
   witness must replay in the full curve group. Require the same ordinary
   outcome and representative counts as baseline. Preserve all failures.
2. Freeze release binaries by hash. Run a dimension-8/10 small screen in
   A/B/B/A order (A = #1421, B = fingerprint candidate), with three
   internal repeats per process. Stop and reject if either dimension's
   two-process median query exceeds 1.20 times baseline. Do not run the
   full panel after a small-screen failure.
3. Only if the small screen passes, run the n83 full-base
   A/B/B/A/A/B/B/A/A/B panel, five processes per arm, with no builds
   between arms and a 120-second process limit. Retain all build/query
   intervals, peak RSS, outputs, stderr and exit codes. Keep the
   candidate only if full-base median target-query time is at least 10%
   lower, maximum RSS is at least 15% lower, and median index-build
   time is no more than 10% higher than baseline.

Index construction is reusable target-independent preparation; query is
target-dependent. This is an exploratory **component** diagnostic on a
contended host, never a controlled complete-call speedup. A complete
n83 higher-arity F6 relation solver, one-target IC online interval and
same-point rho reference are absent, so their costs and ratios remain
unknown. An end-to-end claim would require all of them and an isolated
benchmark receipt under the repository's IC measurement contract.
