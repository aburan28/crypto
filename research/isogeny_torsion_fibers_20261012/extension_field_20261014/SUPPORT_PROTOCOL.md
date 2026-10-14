# Frozen exact support census

Input: `FIXTURES.jsonl`, SHA-256
`1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`.
The source is `examples/isogeny_torsion_extension.rs support`; compile the
committed release example before execution. This finite census does not tune
the already executed paired DLP arms.

For each fixture, enumerate unordered pairs of every signed frozen raw base
point. For each nonidentity sum `S`, count whether `S-aT` lies in the named
prime-order subgroup for `a=0,1,2`; check `[r](S-aT)=O`, record distinct
nonidentity support and whether the frozen target belongs. Separately count
the support of `[h]S` and verify that its degree-three-isogeny image has the
same distinct size. Check the cubic `chi(S+T)=f(chi(S))` including poles for
every enumerated sum. Count distinct `chi`-orbit keys and keys that collide
among different projected subgroup points.

Hypothesis: raw exact support can be empty even where cofactor-projected
support is dense. A torsion shift may add exact raw support, but a shift
cannot create a new projected pair sum because `[h]T=O`. The census tests
that statement on the exact four frozen bases; it does not test the pending
polynomial-system oracle. Any sum, curve, subgroup, or map mismatch fails
closed. No performance conclusion follows from this census alone.
