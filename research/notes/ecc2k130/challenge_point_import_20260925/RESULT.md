# ECC2K-130 public P/Q import into the rotated n=131 model

**Decision: the public point coordinates import directly.** The polynomial
`z^131+z^13+z^2+z+1` and integer bit encoding in the earlier ECC2K-130
challenge implementation are byte-for-byte the same field representation as
the rotated n=131 gate. The independently implemented field and point laws
both check the literal P and Q as on-curve, nonzero, annihilated by the
published prime subgroup order
`q=680564733841876926932320129493409985129`, and satisfying
`tau(R)=[196511074115861092422032515080945363956]R` for each `R=P,Q`.
The curve order recomputed from its F2 trace recurrence is `4q`. No
coordinate conversion is necessary between these two repository models.

The [machine-readable fixture](challenge_points.json) records P, Q,
`[4]Q`, the four actual `Q+T` points, and polynomial plus beta-3 normal
coordinates for each nonzero point. Beta 3 has full normal rank 131.
The independent replay recomposes every normal mask, verifies its trace,
and confirms all four translated points have the same `[4]Q` projection.
The exact public `Tr(x_(Q+T))` values by torsion index `O,(0,1),(1,0),(1,1)`
are **0, 0, 1, 1**. Each torsion point is fixed by Frobenius, so the
correct translated identity is `tau(Q+T)=tau(Q)+T`; the checked
`tau(Q)=[lambda]Q` identity alone must not be substituted for the
translated point.

The source/input manifest was committed at `45ad005` before fixture
generation; its SHA-256 is
`784a0155fe0ebd3401ee195d3f850520db28a47a2148ce082b20811d63cb7488`.
The fixture SHA-256 is
`939fac55c31b7145be425ce09f543feba29aa47846a8355d1f0841f6e594113a`
and the independent [validation](validation.json) SHA-256 is
`d7017da0f121f60a05957ac3dc5b2e43807fab52e84f1db113f68f0851577324`.
`ci_replay.py` regenerated both files byte-for-byte on a fresh temporary
path and passed every independent bit-serial/Fermat check.

This closes a representation prerequisite for an actual-target n=131
PDP exporter. It is not a relation or decomposition of Q; it supplies no
challenge logarithm, rank measurement, solver runtime, or rho comparison.
The next bounded gate should feed all four literal `Q+T` targets to a
complete rational-point S6/S7 exporter, include exceptional/infinity
branches, replay every claimed witness by the full group law, and retain
exhausted/timeout separately. Only then is a measured natural-target
support or cost claim admissible.
