# Actual NIST parameter coverage

The eleven-curve panel is P-256; K-163, K-233, K-283, K-409, K-571;
and B-163, B-233, B-283, B-409, B-571. B-163 is **sect163r2**.
P-256, including the spelling `P-256`, resolves to `p256`.

`icx validate <curve> --json` operates on the **named curve**, including
its exact prime or binary defining polynomial, coefficients, subgroup order,
cofactor and generator. Values are exported as hexadecimal strings and never
narrowed to a machine word. Its result is explicitly a parameter/S3 diagnostic.

## Frozen correctness protocol

The reference inputs are the three source files and SHA-256 hashes recorded
in [nist-parameter-reference.json](nist-parameter-reference.json). The source
revision is `d1114083390ddcfde9a063fa7187767c7e7b92e9`. An independent Sage
implementation checks the field modulus, subgroup-order primality, generator
membership, `[r]G=O`, positive cofactor and the Hasse necessary condition.
This does not independently count each full curve's points.
The native Hasse check uses the same exact squared inequality; it no longer
rounds the square-root bound up and admits impossible boundary cardinalities.

For every curve, the public fixtures are `(u,v)=(1,2)` and `(r-2,r-3)`.
Set `w=(u+v) mod r`, check `[u]G+[v]G=[w]G`, and evaluate the correct prime-
or binary-field S3 on their x-coordinates. Check that the wrong target
`[w+1]G` fails both the exact signed point sum and S3. The second fixture
forces full-width scalar handling. Every target coordinate and all exact
parameters are frozen in the receipt; native CLI tests compare them to Sage.

Each binary curve also has fifteen bounded checks for quadratic point lifting,
both signed generator lifts, the linear and pure-square solver branches,
nonzero constant rejection, Artin-Schreier round trips, and trace-one rejection
(all five NIST degrees are odd). The group-law checks cover both infinity
identities, inverse addition, doubling, scalar zero and `[r-1]G=-G`.
The `x=0` point obtained from the pure-square branch must lie on the curve
and double to infinity. Independent Sage polynomial roots and group operations
provide the reference; native outputs must match both generator y-roots and
all fifteen booleans. These checks cover all ten exact binary parameter sets.

Success requires every check on all eleven curves to pass; any failure is
preserved in the JSON receipt and makes the reference command exit nonzero.
These fixtures are bounded public arithmetic checks. They do not measure
natural relation yield, factor-base membership, independent rank, polynomial
solver performance, linear algebra or logarithm recovery. No speedup, total
attack cost, fitted exponent or ratio to rho is established by this panel.

Reproduce from a checkout with Sage installed:

```sh
sage -python tools/ic_nist_parameter_reference.py \
  --source-revision d1114083390ddcfde9a063fa7187767c7e7b92e9 \
  --output /tmp/nist-parameter-reference.json
cargo test --release --test icx
```

The source hashes, rather than the supplied revision label alone, bind the
inputs. Review regenerated receipts whenever the parameter constructors change.

## Coverage and remaining work

| Panel | Independent Sage parameter/public-S3 checks | Native real-parameter checks | Full-size IC pipeline |
| --- | --- | --- | --- |
| P-256 | Passed | CI gate added; execution pending | Unavailable |
| All five NIST K curves | Passed | CI gate added; execution pending | Unavailable |
| All five NIST B curves | Passed | CI gate added; execution pending | Unavailable |

The existing `icx run` path remains a smaller-analogue experiment. Its result
is always `scaled`, even if a caller raises `--envelope`; increasing an
envelope cannot make the arithmetic backend wider. `--require-named-curve`
fails explicitly rather than substituting an analogue. Family-incompatible
size flags, field degrees above 62 and zero repetitions are rejected before
any run. Named parameters must pass their checks before an analogue is built.

Koblitz analogue selection preserves `a`: K-163 uses `a=1`; the other four
NIST K curves use `a=0`. The earlier `a=1`-first fallback selected a different
member of the family and has been removed. Random binary analogues remain
labelled analogues; their random coefficients do not equal the named B curve.
K-163 defaults to degree 11: its `a=1` degree-13 analogue is not admitted by
the existing subgroup/cofactor rule (`r=79`, `h=106`). Explicit degree-13
requests fail with the preserved coefficient in the error. The `a=0` defaults
remain degree 13. The smoke tests use an admitted analogue of each coefficient.

- [x] Check all eleven exact parameter sets independently with Sage.
- [x] Add native full-width parameter and public point/S3 validation.
- [x] Compare native outputs with frozen independent coordinates in CI.
- [x] Add binary point-lifting, trace and exceptional group-law reference checks.
- [x] Check the Hasse necessary condition without rounding; test both boundaries.
- [x] Preserve Koblitz `a` and prevent misleading analogue run regimes.
- [x] Add an exact-parameter solve guard and reject incompatible size flags.
- [ ] Execute and pass applicable native CI on the final PR head.
- [ ] Implement and qualify a full-width factor-base/decomposition pipeline.
- [ ] Implement and qualify full-width relation coefficients and linear algebra.
- [ ] Establish natural verified relation yield and independent rank on real inputs.
- [ ] Establish complete verified DLP recovery and full-cost comparisons where feasible.

Passing this diagnostic does **not** mark those last four items complete.
No practical full-size IC capability is claimed for the NIST curves.
