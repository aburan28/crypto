# Prime-to-characteristic order panel after the V1 screen

Status: frozen before any new target generation, fiber census, relation run or
timing at the three new orders. This is a versioned correction to the V1
candidate selector; the original V1 screen, output and checker remain intact.

The V1 selector admitted `r=p` at its first three selected cells. That is a
different subgroup regime from the order-229 reference, which has `r != p`.
For a comparable prime-to-characteristic trend, select from the **already
enumerated and independently checked** V1 rows in ascending `p`, requiring
`r != p`, distinct `r`, and `r > 229`. The first three are:

| cell | field | subgroup | `#E(Fp)` | `#E(Fp²)` | cofactor | source/destination j |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| A | 409 | 433 | 387 | 167571 | 387 | 84 / 368 |
| B | 421 | 457 | 387 | 176859 | 387 | 55 / 176 |
| C | 457 | 421 | 495 | 208395 | 495 | 119 / 57 |

The order-229 cell `(p,r,h)=(211,229,195)` is the fixed fourth cell. V1
candidate selection was exploratory, so this V2 choice is not described as a
preregistered result. The *measurements* below are preregistered before any
new panel observation. Rechecking these explicit cells from public formulas
and using a new independent point/fiber census is required before a solve.

## Frozen complete comparison

Use the source `y²=x³−21x+37`, `T=(3,1)`, and Vélu destination
`a'=-9a-30x_T²`, `b'=-27b-70x_T³-42a*x_T` over canonical `Fp2::new(p)`.
Audit every source point, destination image and preimage count; establish
`#E(Fp²)=h r`, exact order `r` of a deterministic first cofactor-projected
generator and its image, and the map identity on all `r²` ordered subgroup
pairs. A separate Rust checker must recompute field and curve arithmetic.

For each admitted cell, freeze 64 distinct source subgroup targets by sorting
the public point-tuple digest with domain
`different-j-prime-to-char-target-v2:13`; keep the scalar list checker-only.
Retain the four relation arms (source projected, projected image, destination
direct, full rational image), direct source BSGS, negation rho and rho A/A.
Use eight signed columns, the same deterministic line and full-image selection
rules as the order-229 run, three arm-rotated rounds, prefixes `N=1,4,16,64`,
and one cold process per arm/round/prefix. Freeze the relation cap as `8r`
trials per target, rho as at most 32 restarts of `20r` steps, and BSGS baby
size `ceil(sqrt(r))`; a cap exit remains a failed or censored child. Each
process has a 300-second wall and 512 MiB RSS ceiling. Do not resample targets
or change a cap after results appear. All 64 scalars in each complete child
must pass independent `[k]G=Q` replay; every relation witness, map, row and
rank must pass the independent checker. Keep raw failures and exact hashes.

Count common construction, setup, relation trials including misses, maps,
linear algebra, recovery, verification, memory and output for every arm.
Report `S_group`, Fp multiplication and inversion axes, and all unpriced work
per target and cold total. Compare with direct BSGS and rho on identical
targets. Run a separate L2-or-better interleaved timing session with A/A and
an independent host performance/replay session; neither is inferred from the
macOS L0 data. At four checked orders, a slope is descriptive with intervals;
an exponent claim needs a complete portable unit, stable domain and a fitted
confidence interval disjoint from the direct reference.

If the full rational census or map check fails, stop before the solve at that
cell and retain its exact failure. If a relation arm caps, preserve the
incomplete record and do not impute a scalar. If all three new cells finish,
the panel may answer the bounded size trend for this family; it cannot by
itself establish a general isogeny speedup.
