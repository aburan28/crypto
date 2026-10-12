# Native GHS structural screener

`ghs_screen` audits an ordinary binary curve

`E: y^2 + x*y = x^3 + a*x^2 + b` over `F_(2^N)`

against every nontrivial factorisation `N = n*l`. It validates the exact
polynomial-basis field with Rabin's irreducibility criterion, calls the shared
`ec_trapdoor::audit_curve` implementation for each factorisation, and emits the
corresponding field and curve tower:

```text
F_2 <= F_(2^l) <= F_(2^N)
C/F_(2^l) -> E/F_(2^N)
```

The report carries the Hess magic number, GHS type, exact descent genus, and
Artin--Schreier compositum degree for every row. Rows are sorted by increasing
genus. The default genus bound of four is only a structural filter.

## Usage

Arguments accept decimal integers or `0x`-prefixed polynomial bitsets. The
modulus includes its leading `z^N` bit; coefficients use bit `i` for `z^i`.

```sh
cargo run --bin ghs_screen -- \
  --degree 8 \
  --modulus 0x11b \
  --a 0 \
  --b 1 \
  --genus-bound 4
```

Use `--output report.json` for a file. Output uses schema `ghs-screen/v1` and
stores genus and cover degree as decimal strings so consumers in languages
without arbitrary-width JSON integers cannot silently truncate them.

For the example above, the tool reports the three factorisations `(n,l) =
(2,4), (4,2), (8,1)`. Invalid coefficient encodings, `b = 0`, reducible field
polynomials, inconsistent degrees, and fields wider than 4096 bits are rejected
before arithmetic begins.

## Exact genus arithmetic

The older `ec_trapdoor` helpers returned `u32` and evaluated
`1u32 << (magic_m - 1)`. That expression cannot represent `magic_m >= 33`.
`ghs_genus`, `ghs_genus_with_type`, `DescentRow::genus`, and this screener now
use `BigUint`; for example, magic numbers 33 and 131 produce the exact generic
genera `2^32` and `2^130`.

## Scope

This tool screens algebraic structure. A row inside the chosen genus bound is
not by itself a vulnerability or speed result. The screen does not compute the
Jacobian order, prove transport of the target subgroup, collect relations,
solve a discrete logarithm, or compare total cost with a matched Pollard-rho
reference. Those steps belong to a separately registered, end-to-end
experiment under this repository's accounting rules.

The reported `C -> E` label describes the GHS descent cover associated with the
field tower. It is separate from the same-field elementary covers cataloged by
[`curve_cover_check`](curves/COVERS.md).

For a checked point-level trace on a subfield-defined binary curve, use
[`ghs_transport`](GHS_TRANSPORT.md). Higher-magic rows remain structural
until a smooth model and norm-conorm map are constructed.
