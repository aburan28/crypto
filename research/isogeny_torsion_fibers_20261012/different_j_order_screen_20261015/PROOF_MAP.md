# Lean arithmetic proof map and geometric boundary

The readable arithmetic argument starts from source coefficients
`a=-21`, `b=37` and the checked rational kernel generator `T=(3,1)`.
Substitution into the degree-three Vélu coefficient formulas gives
`a'=-9(-21)-30·3²=-81` and
`b'=-27·37-70·3³-42(-21)·3=-243`. For each accepted prime,
the native checker reduces these models modulo `p`, evaluates
`j=1728·4a³/(4a³+27b²)` with a nonzero denominator, and obtains the
table in [RESULTS.md](RESULTS.md). The Lean `JCertificate` explicitly
checks a multiplicative inverse for each denominator and the stated
resulting residue. Integer factorization of each direct point count
and the selected batch capacity are also checked.

| Report statement | Lean declaration | Exact scope and status |
| --- | --- | --- |
| The destination coefficient formulas evaluate to `-81` and `-243`. | `OrderScreen.destinationA`, `destinationB` | Proved integer substitutions. The Vélu coefficient formula itself is an input to this arithmetic certificate. |
| Source/destination j residues are `77/160` at 211, `183/222` at 223 and `33/164` at 241. | `p211SourceJ`, `p211DestinationJ`, `p223SourceJ`, `p223DestinationJ`, `p241SourceJ`, `p241DestinationJ` | Proved modular numerator, denominator inverse and result for the six literal model coefficients. No Lean elliptic-curve `j` definition is attached. |
| Each accepted pair has distinct listed j residues. | `p211DifferentJ`, `p223DifferentJ`, `p241DifferentJ` | Proved numeric inequality modulo each named prime. The general theorem that isomorphic nonsingular curves have equal j is not formalized here. |
| `44,655=3·5·13·229`, `50,175=3²·5²·223`, `58,563=3⁵·241`; the selected subgroup has room for 64 nonidentity targets. | `p211OrderFactor`, `p223OrderFactor`, `p241OrderFactor`, `p211DistinctCapacity` | Proved finite arithmetic. The link from the direct native point census to those literal counts is outside Lean. |
| The selected audit counts satisfy `14,885·3=44,655`, `229²=52,441`, and `229 mod 3=1`. | `p211FiberCount`, `p211SubgroupPairCount`, `p211KernelCoprime` | Proved finite arithmetic. The link from pointwise enumeration to these counts is the separate exhaustive Rust runner/checker, not the Lean theorem. |
| `T` is a point of order three; the explicit Vélu map is a homomorphism with exact three-point kernel and the reported j values. | Open | The selected `p=211` instance passed a full two-implementation finite audit of 44,655 points and 52,441 subgroup pairs; the other accepted cells have sampled map witnesses. A faithful Lean `F_(p²)` field, elliptic-curve group law, Vélu construction and general proof remain missing. The finite audit is not a Lean curve-level theorem. |

The project pins Lean `v4.34.1` (toolchain commit
`5045d0056413266e57c625dcd7c365b10e377c52`) with no external
packages. `lake build OrderScreen` from `lean/` passed on the third
attempt; [V3 stdout](LEAN_BUILD_V3.stdout) lists all fifteen
`#print axioms` results as axiom-free. The [first build](LEAN_BUILD_V1.stderr)
used the study root instead of `lean/`; the [second](LEAN_BUILD_V2.stdout)
failed to synthesize `Decidable (JCertificate ...)` before unfolding
the definition. Both failures are preserved. `lean4checker --fresh`
was unavailable. The compiled proof covers the arithmetic rows only.
The separately frozen [full pointwise audit](FIBER_REPORT.md) passed
for `p=211`; the geometric Lean bridge remains open. The three new
finite count equalities are checked in the
[V4 build receipt](LEAN_BUILD_V4.stdout), with their own `#print axioms`
outputs. `lake build OrderScreen` passed, and all eighteen declarations
in V4 report no axioms; [V4 stderr](LEAN_BUILD_V4.stderr) is empty.
