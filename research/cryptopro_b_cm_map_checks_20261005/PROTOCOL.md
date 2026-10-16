# Native correctness replay protocol

Frozen before the native execution on 2026-10-05.

The hypothesis is that an independent Rust replay of the supplied public
CM map reproduces the recorded forward-map correctness outcomes on every
frozen input. This is a correctness experiment; it asserts no performance
gain or discrete-logarithm result.

The input is `endomorphism155.json`, hash-pinned in the verifier. The known
scalars and reference outcomes are frozen in `evidence/legacy/results.json`.
The original Python source and README remain alongside that result. The
Rust replay records both input hashes, uses exact big-integer arithmetic,
and introduces no new random sampling or unknown target.

Success requires all exact component identities, chain and isomorphism
checks, parameter identities, 39 known-input cases, and 16 additivity
pairs to pass, with every recomputed per-input result equal to the frozen
reference. Tests must also reject an altered map coefficient and an altered
final isomorphism. A mismatch terminates the replay with a nonzero exit;
CI preserves the failed job's output. Do not revise frozen inputs to hide
a failure.

The command and complete native report are retained by the dedicated CI
job. No timing is collected. Construction, unsuccessful relation searches,
independent relation collection, linear algebra, and an elliptic baseline
are outside this replay, so total costs and speedup remain unknown. No
scoreboard performance figure or method verdict is changed by this check.
