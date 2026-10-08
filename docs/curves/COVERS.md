# Automatic hyperelliptic cover checks

`curve_cover_check` constructs and replays explicit covers **H -> E** for
the elliptic models in the ICV1 catalog: short Weierstrass models over
`GF(p)` and over extensions `GF(p^k)` with `p > 3`, and ordinary binary
models. Findings are in
[`covers.json`](covers.json), joined to [`registry.json`](registry.json)
by slug **and the full SHA-256 of the exact model JSON**. The registry's
model and EC1 identities stay unchanged. The lab browser shows the finding
and certificate on each curve page.

```sh
cargo test --bin curve_cover_check
cargo run --bin curve_cover_check --
cargo run --bin curve_cover_check -- --check
python3 scripts/build_lab_browser.py
```

`--registry PATH --output PATH` accepts another ICV1 registry. `--check`
recomputes the report, including all equations and checks, and requires
byte-identical saved output. It rejects stale registry/source digests and
altered certificates. No network or external algebra service is used.

## What a finding means

| Status | Meaning |
| --- | --- |
| `verified_over_declared_field` | Explicit polynomial map, nonsingularity conditions, genus and degree checks pass over the specified field. |
| `unsupported` | These constructors do not handle the model or characteristic. Existence stays null. |
| `invalid_input` | Field, coefficient or elliptic nonsingularity checks failed. Existence stays null; the batch retains the row and exits unsuccessfully. |

An identity mismatch, duplicate slug or malformed registry aborts the batch.
A failed candidate never becomes a claim that *no* cover exists. These
constructors are complete for the stated supported model families; they
are not an exhaustive search over every genus, degree or correspondence.
Minimum genus and minimum degree remain null.

Prime moduli pass a fixed 20-base Miller-Rabin screen. This is **not a
primality proof**: that field's primality remains an explicit registry input
assumption, including when a row says `verified_over_declared_field`.
Binary defining polynomials pass the exact Rabin irreducibility criterion,
and so does an extension's modulus `t^k + c_{k-1}t^(k-1) + ... + c_0`, over
`GF(p)`: `t^(p^k) = t`, and `gcd(t^(p^(k/l)) - t, f) = 1` for every prime
`l | k`. An extension's field label must be the one ICV1 derives from `p`,
`k` and the modulus, and its prime `p` passes the same screen as a prime
field's.
Canonical coefficients, field-label consistency and the elliptic
discriminant condition are checked. Field sizes above 4096 bits are rejected.
The checker certifies model geometry, not the catalog's group order,
generator or EC1 subgroup parameters.

The certificate stores ascending-power coefficient arrays for

`H: v^2 + h(u)*v = f(u)` and `x = x(u), y = y_v(u)*v + y_0(u)`.

A field element is one `0x` integer: a prime-field element itself, a binary
element its polynomial-basis bits, and an element
`e_0 + e_1*t + ... + e_{k-1}*t^(k-1)` of `GF(p^k)` the integer
`e_0 + e_1*p + ... + e_{k-1}*p^(k-1)`.

The curves mean their smooth projective models. The verifier substitutes
the supplied map into the target equation, reduces using H's equation,
and checks both the constant and v coefficients. It separately checks
the construction's genus and degree hypotheses; an equation identity
alone would not certify them.

## Prime characteristic greater than 3

For `E: y^2 = F(x) = x^3 + a*x + b`, choose the first
`c in {0,1,2,3}` with `F(c) != 0`. Such a c exists: a cubic has at most
three roots and these four elements are distinct. Construct

`H: v^2 = F(u^2+c)`, with `(x,y)=(u^2+c,v)`.

The verifier checks that the resulting sextic is squarefree, using
polynomial gcd with its derivative. This also follows from E being
nonsingular and `F(c) != 0`: a repeated root would force either `u=0`
or a repeated root of F. Thus H is a geometrically integral hyperelliptic
curve of genus 2. Both H and E are quadratic over their coordinate lines;
the line map has degree 2, so the function-field tower gives
`[K(H):K(E)]=2`. The map is separable and defined over the **same** field.

**Over an extension `GF(p^k)`** the construction is the same, under the
same name (`prime_quadratic_pullback_v1`): nothing in it uses that the field
is prime. The constants `0, 1, 2, 3` are distinct in characteristic
`p > 3`, the cubic still has at most three roots, and squarefreeness is the
same gcd, taken over `GF(p^k)`. The map is defined over the model's own
field `GF(p^k)`, not over `GF(p)`.

This elementary construction is stronger for this one-map question than
the previously discussed genus-at-most-5 result about *two independent*
maps. In particular, prime-field catalog models do admit degree-2 covers.
This does not impose a special CM condition.

## Ordinary binary models

For `E: y^2 + x*y = x^3 + a*x^2 + b`, require `b != 0` and compute
`d = sqrt(b) = b^(2^(m-1))` in its exact polynomial-basis field. Construct

`H: v^2 + u^2*v = u^7 + a*u^4 + d*u`,

with `(x,y)=(u^3,u*v+d)`.

The verifier checks `d^2=b` and the identity

`y^2+x*y-(x^3+a*x^2+b) = u^2 * (v^2+u^2*v-(u^7+a*u^4+d*u))`.

On `u != 0`, put `s=v/u^2`. Then

`s^2+s = u^3 + a + d/u^3`.

The two poles, at zero and infinity, have odd order 3 and nonzero leading
coefficients. They cannot cancel by an Artin-Schreier change of variable;
the quadratic extension is geometrically integral and separable. The
Artin-Schreier genus formula gives `2g+2=(3+1)+(3+1)=8`, hence **g=3**.
The affine model is smooth: the v derivative is `u^2`; at its only zero,
`u=0`, the u derivative is `d != 0`. In the infinity chart
`t=1/u, w=v/u^4`, the equation is
`w^2+t^2*w=t+a*t^4+d*t^7`, which is smooth at infinity.
The coordinate-line map has degree 3; the two quadratic function fields
give `[K(H):K(E)]=3`, also separable in characteristic 2.

This covers the catalog's ordinary binary models, including ECC2K-130.
It corrects any inference that the absence of a useful intermediate
subfield would prevent an ordinary same-field hyperelliptic cover.

## Interpretation and scope

These elementary same-field constructions establish existence. They
perform **no Weil descent to a smaller field**, no subgroup-map computation,
and no relation-solving or discrete-log benchmark. Those fields remain
`not_tested`/null. Existence carries no speedup or vulnerability label.
The elliptic factor already accounts for part of the cover's Jacobian;
its existence alone does not establish cheaper computation there.

The Rust tests enumerate all affine cover points on every nonsingular
short Weierstrass model over F5, F7, F11 and F25 (`t^2 + 2`), on four
models over F125 (`t^3 + t + 1`) and on every ordinary model over F8; every
nonsingular model over F49 (`t^2 + 1`) verifies. Rabin's test over `GF(p)`
counts the monic irreducible polynomials of degrees 2, 3 and 4 over F5,
2 over F7, 4 over F2 and 3 over F3, against the necklace formula. They also reject singular/reducible inputs and altered maps,
coefficients, genus, degree and identity bindings, and replay the entire
catalog deterministically. Point enumeration checks the maps; the genus
and function-field degree arguments above establish the geometric claims.

## References

- Claus Diem, *Families of elliptic curves with genus 2 covers of degree 2*,
  [arXiv:math/0312413](https://arxiv.org/abs/math/0312413), background on
  genus-2 elliptic covers; the elementary pullback proof above is supplied here.
- Arsen Elkin and Rachel Pries, *Ekedahl-Oort strata of hyperelliptic curves
  in characteristic 2*, Algebra & Number Theory 7 (2013), 507-532,
  [doi:10.2140/ant.2013.7.507](https://doi.org/10.2140/ant.2013.7.507),
  Notation 1.1 and the displayed genus formula on page 508.
- Xavier Xarles, *Hyperelliptic curves covering an elliptic curve twice*,
  [arXiv:1303.4220](https://arxiv.org/abs/1303.4220), Corollary 3.
  This broader two-map theorem is context, not the verifier's algorithm.
