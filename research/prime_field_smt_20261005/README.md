# Prime-field decomposition SMT solver

This experiment adds an exact Rust frontend for the finite two-point
decomposition problem on a short-Weierstrass curve over `F_p`:

```text
P1 + P2 = R,  0 <= x(P1),x(P2) < B.
```

It is a decomposition-stage tool, not a complete ECDLP solver. In
particular, a literal prime field has extension degree one, so there are no
extension coordinates on which to perform classical Weil descent. The
three available encodings instead work directly over `F_p`:

- `bit-vector`: portable `QF_BV` with explicitly widened modular
  arithmetic;
- `finite-field`: cvc5 `QF_FF` with full point and slope coordinates;
- `finite-field-s3`: cvc5 `QF_FF` using the Semaev `S3` polynomial and an
  exact trie of liftable factor-base abscissae.

Every SAT result is replayed with independent big-integer curve arithmetic.
For `finite-field-s3`, the checker computes all square-root signs and accepts
only a signed pair whose affine sum is the original target.

## Input

Exact integers are serialized as decimal strings; JSON integers and `0x` or
`0b` strings are also accepted on input.

```json
{
  "schema": "prime-field-smt.instance/v1",
  "curve": { "p": "17", "a": "2", "b": "2" },
  "target": { "x": "13", "y": "10" },
  "x_bound": "4"
}
```

Validation rejects composite moduli, singular curves, noncanonical field
values, off-curve targets, and unsafe exporter widths. For moduli larger
than 64 bits, the receipt labels the fixed-base Miller-Rabin assurance as
probabilistic rather than presenting it as a proof.

## Usage

Build the Rust binary:

```sh
cargo build --bin prime_field_smt
```

Check an input or solve it with the exact native reference:

```sh
target/debug/prime_field_smt check --instance INSTANCE.json
target/debug/prime_field_smt reference \
  --instance INSTANCE.json --max-factor-points 1000000
```

Run the compact Semaev backend with a CoCoA-enabled cvc5 build:

```sh
target/debug/prime_field_smt solve \
  --instance INSTANCE.json \
  --solver /path/to/cvc5 \
  --encoding finite-field-s3 \
  --timeout-ms 120000 \
  --out-dir NEW_RESULT_DIRECTORY
```

Result directories are create-only and retain the frozen input, emitted
query, raw stdout and stderr, hashes, solver arguments, timing, parsed
answer, and verification receipt. `emit` produces a query without running
a solver, and `verify --encoding ...` independently checks retained output.

Native finite-field queries require a cvc5 build configured with CoCoALib;
the official cvc5 1.4.1 non-GPL static binary used in the compatibility
control does not include it.

## Scope

The implementation establishes a reusable exact SAT/SMT backend, but the
frozen 20-bit trial did not complete under the fixed cap. It therefore does
not support a speed claim against Pollard rho, nor an extrapolation to
51- or 83-bit prime-field curves. See [RESULT.md](RESULT.md) for the
recorded outcomes and [PROTOCOL.md](PROTOCOL.md) for the preregistration.
