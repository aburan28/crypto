# P-192 fixed-origin genus-three reduction certificate

**Result.** Every one of the 160 frozen signed Type I quotient classes has a
unique effective degree-three representative relative to the selected origin
divisor (T). The previously checked residual divisor (D_i) has three
distinct geometric abscissae, and its interpolation is
(y=a_i x^2+b_i x+c_i) with (a_i\ne0) in every class. Thus no line contains
the three geometric support points. On the smooth nonhyperelliptic
genus-three plane quartic, Riemann--Roch gives (h^0(D_i)=1), so the
complete linear system (|D_i|) contains only (D_i).

This certifies the fixed-origin degree-three reduction more strongly than
the previous support census: its 160 residual divisors are nonspecial over
the algebraic closure, including the 59 irreducible, 81 mixed-splitting,
and 20 fully split cubics. For this origin and these classes, varying an
effective degree-three representative cannot create new shared-support
relations. A collector must change its classes or origin, use a higher
degree/function-divisor relation, or exploit another geometric identity.

## Certificate and proof map

The [screen](../../nonspecial_support_screen.py) reads the immutable
[160-class support receipt](../../type_i_divisor_support.json), checks all
160 quadratic leading coefficients modulo the P-192 prime, and freezes the
coefficients for the [independent PARI replay](../../verify_nonspecial_support.gp).
The [receipt verifier](../../verify_nonspecial_support_receipt.py) binds the
GP input byte-for-byte to the earlier support record and requires the
exact PARI pass marker. The result is 160 nonzero quadratic coefficients
and zero collinearity witnesses. The earlier Python and PARI support checks
certified all 160 cubic (x)-polynomials squarefree, so their three roots
are distinct over the algebraic closure.

The geometric implication uses the canonical plane-quartic embedding.
A line through all three points would be either vertical, impossible for
three distinct abscissae, or (y=mx+n). In the latter case the quadratic
(a_i x^2+(b_i-m)x+(c_i-n)) would have three distinct roots, forcing
(a_i=0). The [Lean theorem](../../lean/NonspecialSupport.lean)
formalizes this algebraic noncollinearity implication over an arbitrary
field and is fresh-kernel checked with its axioms audited. The smooth
plane-quartic model and concrete coefficients are checked by the prior
PARI/Python receipts. The final Riemann--Roch bridge from noncollinearity
to (h^0(D_i)=1) is a written argument; it is not yet a Lean theorem about
the concrete quartic and its divisor class group.

Specifically, (K) is a line section of degree four. If
(h^0(K-D_i)>0), then (K-D_i) is equivalent to an effective divisor of
degree one, and a line section contains (D_i), contradicting the
certificate. Hence (h^0(K-D_i)=0); with (deg D_i=g=3), Riemann--Roch
gives (h^0(D_i)=1). This uniqueness is relative to the fixed origin
(T). It does not construct a generic Jacobian group law or new target
relations.

## Reproduction

From the repository root:

```sh
python3 research/p192_extension_cover_20261012/nonspecial_support_screen.py research/p192_extension_cover_20261012/type_i_divisor_support.json research/p192_extension_cover_20261012/type_i_nonspecial_support.json research/p192_extension_cover_20261012/type_i_nonspecial_support_input.gp type_i_divisor_support.json
gp -q research/p192_extension_cover_20261012/verify_nonspecial_support.gp
python3 research/p192_extension_cover_20261012/verify_nonspecial_support_receipt.py research/p192_extension_cover_20261012/type_i_divisor_support.json research/p192_extension_cover_20261012/type_i_nonspecial_support.json research/p192_extension_cover_20261012/type_i_nonspecial_support_input.gp research/p192_extension_cover_20261012/type_i_nonspecial_support_gp.log research/p192_extension_cover_20261012/nonspecial_support_screen.py research/p192_extension_cover_20261012/verify_nonspecial_support.gp
```

The pinned Lean project builds `NonspecialSupport`; its fresh checker and
axiom audit are wired into the P-192 proof CI job. The [run manifest](RUN.json)
lists the exact source and receipt hashes and retained commands and statuses.
