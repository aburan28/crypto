# Native ICV1/FB1 reproduction protocol

Date frozen: 2026-10-04

This protocol tests one representation claim only: that the depth-18
Dickson-torus factor base on
`icv1-fp256-t89188191154553853111372247798585809583-f188c491`
(standard name P-256) can be regenerated natively and named under this
repository's ICV1 and FB1 conventions.  It is a storage and correctness
experiment, not an ECDLP speed experiment.

## Hypothesis

A native Rust generator, using `CurveParams::p256()` and the registered curve
identity, will reproduce all of these exploratory targets:

| quantity | frozen target |
|---|---:|
| Dickson depth | 18 |
| root exponent | `0x30000000000000000000` |
| terminal trace | 0 |
| negation-folded columns | 131,239 |
| materialised signed points | 262,478 |
| legacy `SHA256(index || x || min(y,p-y))` | `e52a7b606604641efc99d5a0ef933bb258457e4b65f3c3112be88d792d9e217e` |
| sorted wide point-key SHA-256 | `8fa207da35a4426d3992a1b8e9b08ccbfc8795208f450cf76b44bd47d05bd915` |
| FB1 | `FB1h2ea06bef7f7a` |
| full FB1 SHA-256 | `2ea06bef7f7aa68ce83ad0889116a621885ed03937b5de5481ace9e1dc3626ad` |

The FB1 preimage is the repository's canonical
`ecbench.factor_base/v1` object.  Its builder family is `dickson-torus` and
its parameters are the depth, root exponent, terminal trace, and explicit
point-key encoding
`prime-affine-x-plus-one-shift-sign/be33/v1`.

## Frozen inputs and reference

- repository base commit:
  `2a3b92792c6a78b1ce28d9a89e36fe33e2705ed3`;
- curve identity: the P-256 entry in `docs/curves/registry.json`;
- curve construction: `CurveParams::p256()`;
- norm-one generator: the projection `(1,5)^(p-1)` in
  `F_p[i]/(i^2+1)`;
- two-primary generator: the preceding generator raised to
  `(p+1)/2^96`;
- factor-base fibre: the trace coset rooted at
  `g96^(3*2^(94-18))`.

The target hashes above came from an exploratory, non-native implementation.
They are comparison inputs, not accepted results.  Only the native replay
produced after this protocol is versioned may promote them to repository
evidence.

## Native procedure

1. Recompute the ICV1 identity from the actual P-256 parameters and require
   the registered slug, full ICV1, EC1 alias, and curve UID.
2. Prove the chosen element has exact order `2^96` by checking its
   `2^96` and `2^95` powers.
3. Enumerate all `2^18` trace-fibre abscissae; check the terminal Dickson
   value, lift each quadratic residue, and retain the smaller ordinate.
4. Require distinct sorted abscissae and independently check every stored
   point against the curve equation.
5. Materialise both signs with one shared column and coefficients `1` and
   `n-1`.
6. Hash sorted point keys.  The point key is the repository's prime-field
   rule `((x+1)<<1)|sign`, widened from eight to 33 big-endian bytes;
   `sign=1` iff `y>p-y`.
7. Canonicalise and hash the `ecbench.factor_base/v1` identity object.
8. Write and read back an `ecbench.factor_base_dump/v1-wide` dump, then
   have `ecbench db sql` accept it.  Wide curve integers and coefficients
   are decimal strings because ordinary v1 is intentionally `u64`-limited.

The implementation and every experiment command must be native Rust.  The
full point dump is deterministic derived data; the committed result records
its SHA-256, byte count, regeneration command, and the complete compact FB1
preimage.  A small dump is retained as a parser fixture if needed.

## Success and stop conditions

Success requires every target in the hypothesis table, a passing full dump
round trip, zero off-curve or duplicate points, native unit/integration tests,
and successful database ingestion.  Any mismatch stops the run and is kept as
the result; no target may be changed after seeing native output.

The relation-count corollary may be recomputed exactly for 17 distinct
signed columns.  Its Poisson model and the selector regularity bound remain
mathematical diagnostics.  This protocol measures no decomposition solve,
relation collection, linear algebra, or matched rho reference, so it cannot
produce an end-to-end `S`, an exponent improvement, or a speedup claim.
