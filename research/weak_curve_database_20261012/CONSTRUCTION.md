# Constructing and certifying a covering candidate

The cubic norm-one family gives a direct, reproducible path to an elliptic
curve with the recorded genus-three covering geometry. Work over
L=F_(q³)=F_(p⁶) with q=p² and p>3. Fix and record the irreducible polynomial
for L. Choose a nonidentity element from

    T = {lambda in L× : lambda^(q²+q+1) = 1}.

The forward construction is lambda=alpha^q/alpha with alpha in L×. Its norm
is one because the three conjugate ratios telescope:

    (alpha^q/alpha) (alpha^(q²)/alpha^q)
        (alpha^(q³)/alpha^(q²)) = 1.

For a chosen lambda in T, the finite-field form of Hilbert 90 supplies
alpha satisfying alpha^q=lambda*alpha. Require alpha outside F_q and all
roots in the model below distinct. The elliptic model is

    E_alpha : y² = x (x-alpha) (x-alpha^q).

It has full rational 2-torsion. Record exact coefficients, field encoding,
model hash, both permitted twist signs, certified point count and ordinary
trace. The published covering construction then produces a smooth
genus-three curve C over F_q. Verify its equation, genus, covering map and
induced isogeny to the full Weil restriction Res_(L/F_q)(E_alpha).
Each is a separate certificate; the norm equation establishes the
candidate model.

For odd p, q=p² is one modulo four and the torus order q²+q+1 is
three modulo four. Hence every chosen lambda has the explicit fourth
root lambda^((q²+q+2)/4), checked by the
[Lean fourth-power theorem](lean/WeakCurveDatabase.lean). In the cited cubic branch
proof, a rational 2-isogeny gives a target with full rational 4-torsion.
Its order is divisible by 16, and the trace congruence and
endomorphism-order argument imply the necessary condition v₂(f_pi)>=2
for both allowed twist signs. This screen is useful during class
construction. It is not sufficient: the p41, p43 and p47 diagnostics
retain respectively 396, 376 and 424 ordinary high-depth zero rows
under their randomized census method. The exact two-branch class
decision uses the selected CM-order intersection or the complete
exceptional-unit oracle from the parent ISO-1 report.

To turn the geometric candidate into a selected-subgroup route, specify
the subgroup order r, generator G, target distribution and explicit maps.
Verify that the kernel is disjoint from the chosen subgroup and that
the covering/isogeny map transports G and each target as claimed. Record
CM support preparation, failed searches, model construction, map
construction, map evaluation, field conversion, destination preparation,
and target-dependent work as distinct costs. A comparison uses one
previously unseen public point and a same-point rho baseline under the
repository's measurement contract.

## Recorded worked model

The imported p7 corpus contains the exact model

    ICV1:fpk-7-6-26accc5b:-178:117828:1,2,5,0,0,4:unk:unk:r:26dd618d41b9

over F_(7⁶), with polynomial-basis modulus coefficients
[5,6,5,3,5,5]. Its short Weierstrass coefficients are
a=[0,2,5,6,1,3], b=[1,0,0,2,2,5]. The certified trace is -178, so
the order is 7⁶+1-(-178)=117828. The p7 class oracle marks the cubic
branch supported and the quadratic branch unsupported; the model
index marks this exact cubic model yes, with coordinate replay evidence
linked in curve_model.record_json. Its selected subgroup fields are
unspecified. This is a construction and identity witness. A selected
subgroup comparison requires the additional witnesses above.

## Other construction channels

The machine-readable [mechanisms.json](mechanisms.json) gives eleven
further recipes: anomalous trace-one curves, smooth selected groups,
low embedding-degree pairing transfer, supersingular j=1728 reductions,
the nonsplit quadratic torus branch, Type I and Type II coverings, binary GHS
coefficient control, odd-prime polarization, isogenous covering
representatives, and selected points from a proper subfield. Each names
its recognition test, exact witnesses, subgroup gate and validation
status. Published mechanisms retain their source attribution.
