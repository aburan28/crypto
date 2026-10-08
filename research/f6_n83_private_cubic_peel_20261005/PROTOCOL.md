# n83 F6 five-sum private cubic peel after the quartic certificate

Registered before implementation or measurements. This is an exact
structural follow-on to #1446 at head `58291e8ce`. That result
certified 14,058 of 29,880 one-round degree-four prolonged rows but
left 15,822 unresolved at all 90 source variables. The question here is
whether some survivors can be removed using private degree-three
monomials, preserving every degree-two-or-lower consequence and hence
every source-only affine consequence.

Use the same registered K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, standard dimension-18
factor base (261,447 geometric points, 261,444 subgroup-usable points,
130,722 folded columns), 332-equation five-summand Boolean system,
and public T001 target
`(355fb5df7a905f16921eb,5900a390f42d290f1bbe)` at torsion
offset zero. Use the same 90 source-bit multipliers in bit-major order,
one round only, exact Boolean products, and 32 evenly spaced degree-four
candidates per row. Recompute the parent certificate and require exactly
15,822 unresolved row IDs and frozen BLAKE3 witness digest
`153342de74d4e40a2ea5f96ba35d24b7bdf2726700ed9b04714c9c89ca49dbdd`
before the new cubic scan. Repeat the planted `[0,2,4,6,8]` control,
checking every generated product and full group replay.

For each of the 15,822 quartic-unresolved products, extract its exact
sorted degree-three monomials. Select up to 32 candidates using index
`floor(j * length / count)` as in the parent certificate. Count each
selected exact 512-bit monomial against **all** 15,822 unresolved
products and all 332 original equations in a second pass. Certify a
product row only if one candidate occurs exactly once in this remaining
row set. A uniquely occurring degree-three term forces its row
coefficient to zero in any degree-two-or-lower combination. Original
rows are retained and included in occurrence counts, so the new
unresolved set is conservative. Record the exact certified/unresolved
counts, candidate count, term occurrences, row IDs and ordered witness
digest. The quartic and cubic certificates together preserve every
degree-two-or-lower consequence of the full one-round row span.

If at most 2,000 prolonged rows remain, construct their core together
with all original equations and attempt the deterministic root reducer
with an explicit 6,500,000-column cap. Structural success is a completed
reduction with a contradiction or nonconstant source-only affine row;
then repeat all four torsion offsets and verify any decomposition in the
full curve group before proposing a solver. A completed reduction
without either rejects this two-level peel for this ordinary input.
If more than 2,000 survive, record the count and stop; that result is a
bounded feasibility diagnostic, not evidence of absent constraints.

Run one native Rust release binary with one thread, 7-GiB observed RSS
gate and 300-second process timeout per planted/ordinary process. A
host wrapper samples process RSS and kills above 7 GiB. Preserve raw
stdout/stderr, status and RSS samples, exact code/input/binary hashes,
build and focused-test logs, and a SHA-256 manifest. CPU wall times on
the contended host are exploratory only. This probe does not solve a
natural relation, one target, or a paired rho reference; candidate ID
and F6/IC speedup remain unknown.
