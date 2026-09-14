# Reuse coefficient circuits before trying stronger pruning

This follows the constructor floor in [blocks_01](../blocks_01/BOUNDS_AND_NEXT.md).
The frozen direct constructor cannot cross the matched S3 multiplication count
even with free, perfect pruning. This experiment changes that constructor,
while keeping the exact coefficient order, block rule, leaf solver and verifier.
The [contract](contract.json) is frozen before timing. `RESULTS.md` records
the outcome; algebraic identities below are not speedup claims.

## Why a quadratic table is exact

For a binary subspace V and fixed target abscissa r, the admissible coefficient
set B_r={b : b²+b+r in V} is an affine F2 subspace, or empty: b -> b²+b is
linear. Its size is 2^k when nonempty. Include b=0 while constructing this
affine space; exclude it only during the actual search. Gaussian elimination
finds a basis, and explicit enumeration checks its coverage against the
predecessor's exact Artin--Schreier fibers. No b value is silently discarded.

The flattened support circuit contains C, Lz, Lw and all mixed entries M_ij.
From the predecessor's expansion with h=b²+b+r, the constant and linear
columns depend linearly on b² and b⁴, and

    M_ij = Q_ij + b⁴*(U_ij+r*T_ij)
                   + b²*(U_ij+r*W_ij+r²*T_ij) + b³*r*T_ij,

where T,U,W,Q depend only on the factor-base basis. Frobenius powers b²,b⁴
are binary linear, while b³=b*b² is binary quadratic. Thus the entire circuit
vector is quadratic in the coordinates t_i of b=b0+sum_i t_i*e_i.

Write that vector as C0+sum_i t_i*Di+sum_{i<j}t_i*t_j*Dij. Evaluate it directly
at b0, b0+e_i and b0+e_i+e_j: the unique coefficients are

    C0  = C(b0),
    Di  = C(b0+e_i)+C0,
    Dij = C(b0+e_i+e_j)+C0+Di+Dj.

That requires 1+k+k(k-1)/2 direct circuit evaluations. To materialize the
table, adding coordinate i extends the old half by C(t+e_i)=C(t)+Di+
sum_{j<i}t_j*Dij. Build the latter affine-linear direction table by repeated
doubling and XOR; this needs no new field multiplication. Uniqueness of the
quadratic expansion proves that every table entry equals direct construction.
Small full-table and larger sampled direct comparisons check the implementation.

At k=8, the number of direct evaluations falls from 256 to 37. The direct
constructor uses 259 multiplications at d8, so its sample evaluations cost
9,583 multiplications per target, versus 66,304 before reuse. These counts
exclude affine-space construction, interpolation, materialization, inversions,
rank tests and verification. They are a constructor component, not total cost.

## Two candidates, held apart

`reuse-circuit` uses the table and retains individual inversions of nonzero b.
`reuse-batch-inverse` additionally batches those inversions during table setup.
Both use the unchanged block rule and predecessor b order. All inversions,
table construction and failed work are charged. Batch inversion changes the
cost of inverses, not the solutions or the block certificates.

The table contains 2^k*(d+1)^2 field slots: 20,736 at k=d=8. This is a logical
slot count, not measured resident/peak memory. Python objects, dictionaries,
temporary derivative tables and allocator overhead remain outside that count.
Table materialization uses (2^(k+1)-k-2)*(d+1)^2 field additions, on top of
interpolation and direct evaluation work. A multiplication-only gain can be
cancelled by additions, rank tests or storage.

Precomputing every circuit can worsen time to a first relation. Those
regressions stay in the first-mode comparison; they are not hidden by reporting
complete enumeration alone. The r=0 chart explicitly falls back to the frozen
predecessor.

## Boundaries and experiment

The signed-triple success ceiling is unchanged. The direct S3 pair table is
the strongest implemented reference, rerun on every matched input; chained
S3 and symmetric S4 SAT remain controls. Field add/mul/square counts and selected
binary operations remain separate, with no calibrated total-cost or rho ratio.
The original 20% saving across three increasing sizes is not replaced by a
constructor-only target.

All 16 predecessor inputs and eight fresh holdouts are frozen, yielding 288
three-second cold trials across six variants and both first/enumeration modes.
The identical eight-target n30,d8 predecessor batch is rerun with shared setup
charged and a 60-second whole-batch budget for the reference, both candidates
and S3. Fresh holdouts are in the cold corpus. No further size/budget expansion
is authorized by this experiment contract. It is a family-specific exploratory
screen, not a broad WDSat regression or a complete ECDLP experiment.

Any incorrect table, false rejection, invalid relation or complete-set
disagreement invalidates correctness. Lower constructor multiplication count
alone does not establish scaling or a total-cost win. All source, input and
output hashes, timeouts and regressions must be retained with the scoreboard.

## Reproduction

From the repository root, with Python and pycryptosat available:

    python research/nagao_relations/reuse_01/trial.py --validate-only
    python research/nagao_relations/reuse_01/analyze.py

Timing commands are `trial.py` and `batch_trial.py` without arguments, in a
separate checkout of the frozen source recorded in publication.json. Target
freezing and raw writers refuse overwrites. Evidence is published as lossless,
deterministic gzip; the audit reads compressed JSONL directly.
