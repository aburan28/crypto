# Pre-timing correction

The first candidate build changed the solver recurrence but left the built-in
synthetic collision fixture constructing multiplicative endpoints. `--test`
therefore failed all eight collision-algebra rows while all preceding field,
start-point and symmetry rows passed. No candidate timing was run. Preserve the
first binary, build log and failed test. Before rebuilding, make the fixture
construct additive fixed-point coefficients under `ECC_PACKED_FIXED_ADD`; the
runtime kernel and solver stay unchanged. Regenerate the source manifest and
candidate patch, then rerun the full built-in suite.

The corrected fixture then exposed its own degenerate case on GF(2^13): fixed-add
walks keep the Q coefficient at one, so a planted collision with `c=0, epsilon=1`
has zero solve denominator. Seven curves passed and only that expected algebraic
degeneracy failed. Before a third build, constrain the fixed-add synthetic fixture
to nonzero Frobenius rotations. Preserve the second binary/log/test. Runtime and
solver code remain unchanged.
