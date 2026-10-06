# n83 F6 five-summand dimension-18 direct-system gate

Registered before the candidate probe and measurements. The preceding
dimension-16 six-summand direct S3 chain admitted an exact 428-variable
system, but own-degree reduction gave no source-only affine condition on
public T001. This gate tests a different arity/base trade: five source
points from the standard dimension-18 subspace. The expected geometric
source count is about 2^18; the **nominal** five-multiset capacity
`C(2^18+4,5)/r = 4.267` uses subgroup order
`r = 2417851639230796216685689`. This is a counting heuristic, not an
observed relation yield or an actual factor-base count.

Freeze the registered K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, its exact field polynomial
`z^83+z^45+z^2+z+1`, and public subgroup target T001
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)`.
Use the existing 512-bit native Boolean system and exactly four S3
links. The five planted source points `[0,2,4,6,8]` are a correctness
control: verify all equations, exact full-group sum, cofactor-four
projection, and the unique rational 4-torsion bridge. This control is
not ordinary relation yield.

For the ordinary public T001, construct the four source preimages
`[4^{-1} mod r]T001 + T_i`, where the `T_i` are the four checked rational
4-torsion points in the existing six-summand probe. For each offset,
record source-system variables, equation and monomial counts, an
own-degree root reduction under the existing 1,500,000-column cap,
rank, contradictions, all affine rows and their support, XOR count,
resource observations, raw stdout/stderr, and exit code. Use one
process per offset with a 120-second timeout and observed 7-GiB RSS
guard; an unsupported kernel memory cap is not silently reported as
enforced. Freeze source/input hashes and all statuses.

The structural admission gate passes only if planted equations and
group/cofactor replay are exact and all four ordinary systems construct
within the limits. A search-pruning gate passes only if at least one
ordinary offset yields a contradiction or a **source-only** affine row
after the own-degree reduction. A column-limit or timeout is
inconclusive. If structural admission passes but pruning fails, archive
the negative result and do not extend this root-only search strategy.
Do not infer natural relation yield from a planted control or a mere
affine row involving free intermediate bits.

This is a bounded algebraic diagnostic. No candidate ID, F4/F5/F6
complete-call or IC/rho speedup, or one-target online wall-time ratio can
be claimed without an ordinary verified relation and complete solver.
This unisolated host also cannot promote CPU wall-time ratios.
