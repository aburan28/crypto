# Follow-up: external variables do not measure native cost

Observed in the exact interface control:265 declared coordinates became396 native internal/outer coordinates, an increase of131, with zero auxiliary decisions. All128 expected models and every declared functional coordinate were verified. This is not a timing comparison.

Primary source: the frozen modified core function Solver::add_every_combination_xor, vendor/cms-primary-cdcl-fork-v2/source/src/solver.cpp, retains the original XOR for Gaussian use while also constructing cut-CNF clauses and fresh internal coordinates. The exact post-cleaning full norm expansion is unmeasured. The control initialized no Gaussian matrix.

New falsifiable variant: cap affine cone definitions by their post-normalization support size and retain shared linear DAG intermediates above the cap. Compare limits8/16/32 and the original DAG against dense cones on the same literal predicate. Charge all declared and internal coordinates, rows, coefficients, clauses and replay work. Preserve all complete tiny primary model sets before hard use. Preregister any later hard selection on new inputs; do not select a performance winner on the same opened holdout.

Exact output contract: every retained intermediate has its own exact parity definition, so eliminating or retaining a linear intermediate is a bijective functional extension. AND gates remain nonlinear exact clauses. Full primary assignments and their original literal norm replay remain mandatory. No original S3 all-roots or end-to-end gain follows from the representation.

Next bounded experiment: a90-second,4GiB,one-core representation/control profile. The existing60-second hard timeout is not retried unchanged.
