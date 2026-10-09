# Characteristic-two composite-field route

The repository now exposes two separate stages for an ordinary binary curve
`E/F_(2^N): y² + xy = x³ + ax² + b`:

1. `ghs_screen` validates the field and enumerates every tower
   `F_2 ⊂ F_(2^l) ⊂ F_(2^N)` with `N = n*l`. It reports the Hess magic
   number, GHS type, genus, and Artin–Schreier cover degree. See
   [GHS_SCREEN.md](GHS_SCREEN.md).
2. `ghs_transport` executes the **trace branch** when both `a` and `b`
   are fixed by the tower Frobenius. It checks the supplied points and
   subgroup annihilator, computes `Tr(P)` and `Tr(Q)`, and reports whether
   the generator survives. This is the executable characteristic-two
   construction currently available for composite extension fields.

For a sample field `F_(2^6)`, the modulus `z^6+z+1` is irreducible.
The point `(0,1)` on `y²+xy=x³+1` has order two; the degree-three trace
to `F_(2²)` keeps it nonzero. The command below is a transport fixture,
not a cryptographic-size subgroup:

```sh
cargo run --bin ghs_transport -- \
  --degree 6 --modulus 0x43 --a 0 --b 1 \
  --relative-degree 3 --p-x 0 --p-y 1 \
  --q-x 0 --q-y 1 --order 2
```

The JSON schema `ghs-transport/v1` gives the image coordinates and one of
`nonzero_subgroup_image`, `generator_killed_by_trace`, or
`curve_not_defined_over_subfield`. The command does not certify that the
claimed annihilator is prime or that `Q=[d]P`; it marks both facts as
unverified. If the order is independently certified prime, `P` is a
nonidentity point killed by that order, and its trace image is nonzero,
the trace is injective on `⟨P⟩` and preserves the unknown logarithm.

The new trace check requires **both** coefficients to be in `F_(2^l)`.
Earlier `descend_m1` checked only `b`; that was insufficient because
Frobenius need not preserve the curve when `a` lies outside the subfield.
The symbolic `m=2` constructor now also requires subfield-fixed `a`, and
the previously supplied affine formula requires its stated `a=1` case.

## 113- and 192-bit field assessments

`sect113r1` uses `F_(2^113)`, a prime absolute extension degree. Its
only proper subfield is `F_2`, and the exact `ghs_screen` run with the
SEC 2 parameters below found a magic number of 113, a genus of
`2^112 - 1`, and no genus-4 candidate:

```sh
cargo run --bin ghs_screen -- \
  --degree 113 --modulus 0x20000000000000000000000000201 \
  --a 0x3088250ca6e7c7fe649ce85820f7 \
  --b 0xe8bee4d3e2260744188be0e9c723 \
  --genus-bound 4
```

This is a structural result for [SEC 2's sect113r1 parameters](https://www.secg.org/SEC2-Ver-1.0.pdf),
not a complete ECDLP cost measurement. [Hess's GHS analysis](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/315278965D8D277A85812D9498236A88/S146115700000108Xa.pdf/generalising-the-ghs-attack-on-the-elliptic-curve-discrete-logarithm-problem.pdf)
also distinguishes composite from prime extension degree.

The transport command accepted SEC 2's published `sect113r1` base point
and subgroup order, verified the point and `[order]G=O`, and returned
`curve_not_defined_over_subfield` for the only tower `N=113, n=113, l=1`:

```sh
target/debug/ghs_transport \
  --degree 113 --modulus 0x20000000000000000000000000201 \
  --a 0x3088250ca6e7c7fe649ce85820f7 \
  --b 0xe8bee4d3e2260744188be0e9c723 \
  --relative-degree 113 \
  --p-x 0x009d73616f35f4ab1407d73562c10f \
  --p-y 0x00a52830277958ee84d1315ed31886 \
  --q-x 0x009d73616f35f4ab1407d73562c10f \
  --q-y 0x00a52830277958ee84d1315ed31886 \
  --order 0x0100000000000000d9ccec8a39e56f
```

`Q=G` here is an input-validation fixture. The command does not infer
the order's primality from an annihilation check; SEC 2 specifies that
parameter separately.

The binary field `F_(2^192)` has proper subfields and can be screened.
The polynomial `z^192+z^7+z²+z+1` is irreducible; `ghs_screen` accepts it
and enumerated 13 towers for the curve `a=0,b=1`. All 13 rows had magic
number one, so the genus value alone says nothing about the order of a
transported subgroup. The exact screening command is:

```sh
cargo run --bin ghs_screen -- \
  --degree 192 \
  --modulus 0x1000000000000000000000000000000000000000000000087 \
  --a 0 --b 1 --genus-bound 4
```

A trace transport can be executed on a
subfield-defined curve if a checked subgroup point is supplied. That
does **not** turn NIST P-192 into a binary curve: P-192 is over a prime
field, and the odd-characteristic ISO-1 cover remains a separate path.
The `ghs_transport` command above, with `degree=192`, `relative-degree=3`,
and the same order-two point `(0,1)`, returned `nonzero_subgroup_image`
over `F_(2^64)`. This verifies 192-bit **field arithmetic and transport**;
it does not exercise a 192-bit prime-order subgroup.

## Remaining construction gates

For magic number at least two, the binary GHS code constructs structural
Artin–Schreier data and an affine model in one case. It does not yet produce
the smooth higher-genus curve, the norm-conorm map on points, or a
large-scale Jacobian index-calculus solve. The trace command refuses to
represent these as an executable transport. The odd-prime Joux–Vitse
implementation supports its specific degree-six tower `F_(p^6)`; no
general odd-characteristic cover for every composite degree is established
here. The arithmetic extension in
[WIDE_ORDER_113_192_ASSESSMENT.md](isogeny-walk/WIDE_ORDER_113_192_ASSESSMENT.md)
does not change those construction requirements.

| Requested capability | Current evidence and state |
| --- | --- |
| Characteristic-two composite-field structural screen | Available through `ghs_screen`; exact 113- and 192-bit field runs above |
| Characteristic-two point transport with wide field elements | Available for subfield-defined trace cases; checked API and `ghs_transport` CLI, with degree-6 and degree-192 runs |
| `sect113r1` | Its published point and annihilator pass validation; only tower over `F_2` has large genus, and trace conditions fail |
| Higher-magic binary GHS cover and index calculus | Unresolved: smooth model, norm-conorm point map, subgroup preservation, solver and costs |
| Odd-prime composite degrees other than six | Unresolved: a cover and transfer construction for those degrees is absent |

The raw JSON from the four 113- and 192-bit field commands is retained in
[`ghs-transport-evidence/`](ghs-transport-evidence/). These are deterministic
field and transport receipts. The order-two 192-bit-field fixture and the
`Q=G` standardized-curve check do not measure attack work on an unknown
logarithm.
