# n83 F6 five-sum degree-four private-monomial core

Registered before implementation or candidate measurements. This is an
exact structural follow-on to #1435 at parent head
`2745ee13a996a26a2dbb905c18ac90539ee88152`. That experiment reached
1,463,392 columns with four source multipliers and a 1,500,000-column
cap with eight, without a new source-only affine row. This gate avoids
materializing the full degree-four column set by asking which prolonged
rows can *possibly* participate in cancellation of all degree-four terms.

## Frozen mathematics and inputs

Use registered K0 curve `icv1-f2m83-tm6151469093347-debefd74`, the
same exact subgroup, the standard dimension-18 source subspace, and the
same public T001 point
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)` with torsion
offset zero. The source has 261,447 geometric points, 261,444 distinct
subgroup-usable points and 130,722 sign-folded columns. The exact
five-summand system has 339 Boolean variables, 90 source bits and 332
original cubic equations. Use the planted source indices
`[0,2,4,6,8]` as a correctness control; every generated product must
vanish under its full coordinate assignment and the group sum must
replay. Planted controls do not estimate natural relation yield.

The 90 source variables are ordered by bit position, then summand:
`v = bit + 18*summand` for `bit=0..17`, `summand=0..4`. The first 16
match #1435. Execute `k=16`, then `k=90` if the first process passes
its resource gate. For each `k`, form the exact Boolean product of
every **original** equation with each selected source variable, using
`x²=x` and XOR cancellation. The original equations remain degree at
most three; newly added rows are not multiplied again. Row ID is
`variable_position * 332 + original_equation_index`.

## Exact private-column certificate

For each prolonged row, take its exact sorted list of degree-four
monomials. Select up to 32 candidates at indices
`floor(j * length / count)` for `j=0..count-1`, where
`count=min(32,length)`. A row with no degree-four monomial has no
candidate. Use exact 512-bit monomial keys, not a truncated hash.
Pass one records the selected candidates and exact row/term counts.
Pass two regenerates **every** prolonged row and counts occurrences
of every selected candidate across all rows. A row is certified only
when at least one of its candidates occurs in exactly one prolonged
row. Such a private degree-four term forces that row's coefficient to
zero in any linear combination of all original and prolonged rows whose
degree-four part vanishes. Therefore certified rows cannot contribute
to a degree-three-or-lower consequence; removing them preserves all
possible source-only affine consequences. Record a deterministic hash
of row IDs and witness monomials and a list of unresolved row IDs.

Run a small synthetic test with unique and shared degree-four terms.
For both n83 `k` values, report total and certified prolonged rows,
rows lacking degree-four terms, candidate keys, degree-four term
occurrences, unresolved row IDs and peak RSS. A 7-GiB observed RSS
gate and 300-second process limit apply; retain any timeout/OOM/cap
as an inconclusive row. Do not infer absence of source constraints
from a partial certificate.

If the full `k=90` run certifies **every** prolonged row, this exactly
proves that one round of source-variable degree-four prolongation adds
no degree-three-or-lower consequence beyond the original system. If
at most 2,000 rows remain unresolved, build a core from the original
equations plus only those unresolved products and run the existing
deterministic root reducer under its 1,500,000-column and 300-second
limits. The core reduction is an exact preservation of all possible
low-degree consequences. A contradiction or nonconstant affine row
involving only source bits is the preregistered positive structural
outcome; if found, repeat all four torsion offsets before proposing
a solver. If over 2,000 rows remain, record the core size and stop.

Build one native Rust release binary and freeze its hash before running.
Preserve source/input/binary hashes, raw stdout/stderr, exit codes,
memory and timing. CPU times on this Mac are exploratory only. This is
not a complete F6 decomposition, IC candidate or IC/rho comparison;
the complete one-target online interval and same-point rho reference
remain unknown until both methods solve and independently verify the
same public target under the repository measurement contract.
