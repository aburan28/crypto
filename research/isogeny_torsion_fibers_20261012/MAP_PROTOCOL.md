# Explicit-map control (registered after census, before this run)

Date: 2026-10-12. Frozen input: `p=13`, source `E: y²=x³+2`, quotient
`E': y²=x³+11`, subgroup `E(F₁₃)` of order 19. Both exact orders were
recorded by the census. The order-3 geometric kernel is at `x=0` and is
Frobenius-stable, but has no nonzero rational points.

Check the standard 3-isogeny formula, for `x ≠ 0`,

```
X = (x³ + 4b) / x²
Y = y (x³ - 8b) / x³,      b=2.
```

The source identity maps to the target identity; kernel points at `x=0`
map to the target identity over the algebraic closure. For every rational
source point, require an image on `E'`, require all rational images to be
distinct, and check `φ(P+Q)=φ(P)+φ(Q)` for every ordered pair of rational
points. Recount both curves directly. Fail on any mismatch. Stop after
this one named public instance; no timing or target logarithm is run.

The formula and bounded checks establish an executable map on this exact
instance. The algebraic identity over arbitrary fields and extension-field
kernel behavior remain supported by the written derivation, not by this
finite check. Complete DLP cost and comparison to matched rho remain open.
