# Planned exact alternative to unsupported large-characteristic resolutions

Status: a reasoning-tier plan, not an executed regularity experiment or an
established all-case result. The stopped pilot and card assessment remain
immutable under sessions/regularity_probe_1820 and cards/. Singular std works
over (integer,p), but mres and regIdeal require coefficient-field capabilities
that this backend does not supply. A native certificate route can avoid them.

## The invariant and the checks

Let S=F_p[x,y,z], with p=72057980207847001. Let A,B be the degree-12
homogenizations of the two split factor-base polynomials, and f the degree-4
homogenization of the supplied f3 polynomial. J=(A,B,f) means the ideal of
these separately homogenized generators, not saturation or homogenization of
the whole affine ideal. Report reg(J) and reg(S/J) separately.

Use only the public mathematical systems in validation/singular_all.sing,
data/ic_factorbases.json and data/family.json. Its 14720 system lines are
ordered curve, generator, then two controls, twenty setup, ten target. Parse
only F,G,H before I=F,G,H, ignoring the returned basis E. Do not read private
IC journals or sampling coefficients. Recompute the factor polynomials,
assert twelve distinct roots, verify the supplied f has degree four with
unique top term x^2*y^2, and enumerate its full 12x12 grid of zeros.
This independently checks the exact hypotheses below on every system.

A and B have no common point at infinity: at z=0 they are x^12 and y^12.
Their projective intersection is the reduced 144-point affine grid. The
three generators of J are minimal: f has degree four; A and B have distinct
leading pure powers not reducible by leading term x^2*y^2 of f.

## Empty grid-zero set

If f has no grid zero, J has height three and is a regular sequence of
degrees 12,12,4 in the three-variable polynomial ring. Its Koszul resolution
gives reg(S/J)=25 and reg(J)=26. This is a certificate-derived invariant,
not a newly observed GB solving degree or a claim that affine (1) has that
regularity. Every random setup/target must pass the empty-grid check; do not
assign a value to an unsupported zero pattern.

## One removed point: intended double-point control

For the sole zero (a,a), let X=x-a*z, Y=y-a*z, U=A/X, V=B/Y and T=U*V.
The residual-grid ideal K=(A,B,T) is the ideal of maximal minors of columns
(-V,0,X) and (0,-U,Y). Its Hilbert-Burch shifts are generators 12,12,22
and syzygies 23,23. Check the three minors exactly and the fixed homogeneous
degrees. Check f*T=q*A+r*B by native polynomial division and record the
quotients with degree at most 14. This proves K=(A,B):f: K is the saturated
reduced residual-grid ideal (height-two perfect, correct finite grid zeros).
Mapping the multiplication-by-f short exact sequence gives the quotient
resolution shifts:

- homological degree 0: [0]
- 1: [4,12,12]
- 2: [16,16,24,26]
- 3: [27,27]

No adjacent shifts agree, so the mapping cone is minimal. Thus the proposed
certificate value is reg(S/J)=24 and reg(J)=25, conditional on all checks.

## Two exchanged zeros: intended distinct-point control

For exactly (a,b),(b,a), a!=b, put W(t)=product over the other ten roots,
U=W(x,z), V=W(y,z), X=x-a*z, Y=y-a*z, s=a+b, and T=(x-y)*U*V.
Use Hilbert-Burch columns

(-V, U, x+y-s*z) and (Y*V, -X*U, -X*Y).

Their maximal minors are T, -B, A (with row-pair ordering 12,13,23).
Check those identities and degrees exactly over the actual prime field.
K=(A,B,T) has generator shifts 12,12,21 and syzygies 22,23. Check f*T is
in (A,B), retaining homogeneous degree-13 quotients. The same residual-grid
and colon argument applies, and the minimal mapping-cone shifts are

- 0: [0]
- 1: [4,12,12]
- 2: [16,16,24,25]
- 3: [26,27]

The conditional certificate values are again reg(S/J)=24 and reg(J)=25.
An extra grid zero invalidates this case classification; retain failure.

## Independent validation and execution still required

Before applying these values, implement the native C++ parser/certifier,
record all hypotheses, actual arithmetic counters and CPU/wall/resource
boundaries, freeze a new protocol/binary/inputs, and use the director/executor
card workflow. Retain all 14720 per-system outcomes and 920 explicit control
membership/minor certificates, plus setup per curve/generator charged once.
Verify the coefficient identities independently, e.g. PARI exact polynomial
arithmetic. No source or result from the stopped pilot is overwritten.

A separate small-prime regression can test the resolution-shift deductions
with Singular's genuine coefficient-field backend, explicitly without
substituting its outputs for actual-curve measurements: p=1000003, roots
0..11 for both factor polynomials, and affine forms x^2*y^2+x^2+y^2+1
(empty zero set), x^2*y^2+x^2+y^2 (sole zero (0,0)), and
x^2*y^2+(x+y-1)^2 (zeros (0,1),(1,0)). On that integer grid each nonzero
value is below p. Test actual minimal resolution shifts against the above
formulas and convention control (x^2,y^2,z^2), reg ideal=4. This regression
is planned, not executed. It needs its own frozen command card.

Even a completed regularity certificate phase does not recover any random
ECDLP scalar or resolve the separate IC relation-yield obstruction.
