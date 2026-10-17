# Direct Semaev root and quotient controls: frozen protocol

Frozen before implementation and execution on 2026-10-15. This extends
[Gate A](../NEXT_GATE_20261015.md) without changing its four public cells.
The input is `../extension_field_20261014/FIXTURES.jsonl`, SHA-256
`1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`;
the earlier χ run specification has SHA-256
`fe7ff00ad80da0e391e0d4d0602b51806acce925e546959771d8d6074480c1e2`.
Keep `p=61,271,787,3319`, subgroup orders `19,61,211,811`, all four
public generators and targets, signed projected bases, three seeds per
cell, line construction, trial caps `min(20r,10^6)`, memory cap 512 MiB,
and independent scalar replay. Answer scalars are read by the checker,
never by either root oracle.

## Hypothesis and fixed reference

For each known signed base abscissa, solve the actual cubic-sum equation
`S3(x_i,x_j,x_R)=0` as a quadratic in `x_j` over the frozen `F_(p²)`.
Root lookup and exact group addition may yield the same relation rows as
projected-pair indexing at a different complete cost. The quotient arm
applies the degree-three map to the base and query first, then solves the
same direct polynomial on destination abscissae. This tests whether the
quotient changes root support after all map and correction work is charged.
The same-target direct negation rho and six-class rho controls, raw-pair,
projected-pair, χ orbit, χ finite-scan, source isogeny projection, and
scaling isomorphism remain in the paired comparison. The derived reference
for a complete method is the cheaper applicable direct method in one
counted operation unit; local CPU equivalents are descriptive until that
conversion is available.

## Exactness before relation timing

For `p=61,r=19`, enumerate the 19 subgroup points, infinity, and the
two torsion exceptional points. For each target, compare the *set* of
unordered signed-base index pairs from each new oracle against independent
signed-pair group addition. Include repeated summands, both signs,
vertical chords, and quotient-map infinity. Check the root polynomial
by substitution, every returned pair by exact source addition, every
mapped point by destination membership, and the quotient homomorphism on
the tested pairs. Any missing or false pair stops the timing run.

## Full run and accounting

Run both new arms for all four cells and three rounds in the existing
pipeline, with identical walk seeds and trial limits. Charge every
polynomial coefficient evaluation, discriminant, square-root attempt,
inversion, failed root, lookup, map evaluation, candidate addition,
factor-base setup, relation check, rank update, linear solve, target
descent, scalar verification, failed trial, CPU, wall time, and peak RSS.
The 10^7 polynomial-call cap applies per arm and cell; mark a hit as
censored. Retain raw JSONL and stderr without replacing older results.
The independent checker must cover every arm/round slot and replay each
recovered scalar with separate field and curve arithmetic.

Stop after three rounds per cell and one independent replay. Compare
relation count, useful rank, query trials, and whole current-thread CPU
with the paired controls. Identical relations plus higher total cost
reject this implementation on these four cells. A local win remains a
bounded measurement and requires portable counted `S`, isolated
interleaving, and second-host replay before any frontier update.
