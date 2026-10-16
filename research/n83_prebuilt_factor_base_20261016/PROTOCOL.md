# N83 retained factor-base review and construction measurement

Frozen before the new run on 2026-10-16. The question is whether a factor
base for degree 83 can be built ahead of the relation queries, retained with
an exact identity, and reused; and what the retained sizes and symmetries
imply for two-, three-, and four-summand full-smooth relations.

## Inputs and comparison

The primary curve is `icv1-f2m83-tm6151469093347-debefd74`,
`y^2+xy=x^3+1` over the polynomial basis with low modulus terms
`[0,1,2,45]`, subgroup order `2417851639230796216685689`, cofactor 4.
The `a=1` curve is a separate diagnostic with its own subgroup order.
Use the native `koblitz_n83_factor_base_export` example without source
modification. Construct the full declared 54-object panel: two curves,
sequential/hash/Gray-prefix public x policies, seeds 2026100801..03, and
signed-Frobenius orbit-column counts 64, 256, 600. Each object has 166K
points. The archived 2026-10-08 panel and its generic replay are the
historical reference; compare exact point-set hashes for matched recipes.

Run three fresh panel constructions, each with a 300-second internal cap, in
separate directories. Preserve every object and manifest if a run stops.
Run the native `replay` on each complete panel in a separate process. Record
outer wall, user/system CPU, peak resident memory where available, process
status, manifest elapsed time, each object construction interval and byte
count, and panel and source hashes. Read-only replay costs include compressed
object read, hashing, decoding, and generic point/label checks. They are
reported as a conservative verified import path, not as an equal-work speed
comparison with construction. Keep build time, import/verification time,
index construction, rank, relation search, and target online time distinct.

For m = 2, 3, 4, independently recompute the exact multiset count
`binomial(166K+m-1,m)` and its uniform nonidentity-target Markov ceiling
`min(1,M/(r-1))` from the frozen panel values. Cross-check the historical
`pilot-01/support-moments.json`. This is a structural support bound only;
do not estimate actual yield from it. Record orbit-column count K separately
from point count, and compare the historical structured dimension-12 base and
the distinct `a=1` four-summand cold work without mixing curve identities.

## Admission and stopping

The new run is a local L0 construction/import measurement, with no claim to
an isolated wall-time speedup. A retained object is admitted only if its
manifest and native replay pass, the point-set hash matches the archived
recipe, and the curve, subgroup, orbit phase and labels match. An incomplete
or capped panel is a preserved bounded row. A failed replay is a correctness
failure. The three repetitions are descriptive; report their spread.

No m = 2, 3, or 4 relation throughput, useful rank, full pipeline cost, or
matched rho result follows from this run. Advancing a selected base for any
of those arities requires a frozen same-instance solver comparison with
natural queries, independent relation and rank replay, full setup and online
accounting, A/A controls, and the repository's isolation gates.
