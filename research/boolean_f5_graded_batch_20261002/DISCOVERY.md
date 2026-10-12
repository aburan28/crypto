# Qualified exact graded Boolean F5 batch discovery

**Decision: reject the registered >2× complete-F5-batch gate.** This was a
complete, resource-qualified discovery on the frozen source. The untouched
holdout seeds were not run. The candidate improved every primary n=24 batch,
but none of the four primary groups reached the preregistered paired median
and 95% bootstrap lower bound above 2.0. The eight n=16/20 batch-32
nonregression groups all passed their lower-bound-above-0.95 gate.

| n=24, batch 32 | Seed | Paired median baseline/candidate | 95% lower | A/A 97.5% floor | Counted GF(2) word-XOR ratio | Gate |
|:--|--:|--:|--:|--:|--:|:--|
| Independent affine | 20261020 | 1.838 | 1.832 | 1.011 | 5.622 | Reject |
| Walk affine | 20261020 | 1.797 | 1.768 | 1.010 | 5.650 | Reject |
| Independent affine | 3141691 | 1.829 | 1.807 | 1.012 | 5.636 | Reject |
| Walk affine | 3141691 | 1.784 | 1.730 | 1.016 | 5.728 | Reject |

The [sealed bundle](qualified_discovery_01/manifest.json) is the canonical
record. It contains [all 48 cells and 18,816 timed calls](qualified_discovery_01/raw.jsonl),
the [native verifier's result](qualified_discovery_01/results.json), the exact
binary/source/protocol, build and test logs, the CPU manifest, the readiness
samples, and the isolation receipt. All manifest member SHA-256 hashes were
read back from the downloaded GitHub Actions artifact and passed. The native
bundle replay in the workflow passed. The run was
[GitHub Actions 37045285547](https://github.com/aburan28/crypto/actions/runs/37045285547)
at commit `53f26c669852108c8128c16341f722a1c49e4f5d` on Linux x86-64,
AMD EPYC 9V74, with logical CPUs 2 and 3 reserved as one physical core and
`RAYON_NUM_THREADS=1`. The worker wall time was 976.835 seconds; its receipt
was uncontended, with 7.05 other-process CPU seconds, preflight PSI some
avg10 of 3.46, and maximum RSS of 368,660 KiB. The readiness gate waited
through seven initially high-PSI samples before accepting its eighth sample.
The raw SHA-256 is
`a4ed20fec5066e8d80d0e187c23275bc6f7e3dc61ffb81ad10ba7e9a8727a7c3`;
the manifest SHA-256 is
`651b2ae2932bcebd2bf14d012c41477671c3b69e9e320393762a90a2eececbe5`.

Across the complete campaign, 4,151 candidate calls used the cache and 553
took the charged fresh-F5 fallback. In each n=24 independent-affine primary
batch there was one fallback; the walk-affine batches had two. The four n=24
groups averaged about 2.00–2.06 seconds for a fresh baseline batch and
1.11–1.14 seconds for the candidate. Cache compilation averaged about
70–72 ms per cold candidate batch, with returned-object destruction around
27–30 ms. These are descriptive means from the raw receipts; the table's
paired medians and bootstrap bounds are the registered decision statistics.
The roughly 5.6× reduction in counted GF(2) word XORs did not translate into
2× complete-batch wall speed because output construction, fixed setup and
fallbacks remain charged.

The mathematical next test is a **support-aware continuation** with a fresh
protocol and fresh seeds. When every column of the word-aligned high prefix
is present, its elimination transform is unchanged even if lower-degree
columns are absent. For a sparse-support assignment, apply the retained
transform to its low suffix, project onto the exact occupied lower-column
sequence used by the source's sorted-column fallback, and resume M4RI at the
same high-prefix word boundary. This must be checked against the exact
ordered reference output on every new fixture. It can remove the current
base and sparse-walk fallbacks; a separate unpack/compilation improvement
may still be needed to clear the wall-time gate. No source tuning or seed
selection is applied to this rejected discovery, and the reserved holdout
remains unused.

This is a finite public Boolean matrix-F5 engineering result. It does not
measure natural Semaev relation yield, independent relation rank, full index
calculus cost, or an automorphism-discounted Pollard-rho comparison. Those
fields remain null in the native result.
