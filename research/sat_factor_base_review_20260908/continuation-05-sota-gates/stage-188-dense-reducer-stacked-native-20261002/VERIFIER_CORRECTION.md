# Stage 188 verifier correction

The first native screen completed both backend children and passed every
terminal and structural correctness check.  The runner then returned exit 2
because its independent assessment replay compared recomputed floating-point
ratios with Rust's exact derived `PartialEq`.  The stored and replayed ratios
can differ at serialization-scale rounding even though they derive from the
same authenticated wall, CPU, and RSS measurements.

This is a verifier defect, not a candidate or measurement result.  The original
`development/screen/verification.json` is retained with its
`paired assessment replay` failure.  The additive correction compares only
ratio fields with relative tolerance `1e-12`; all discrete fields, pair order,
repeat identities, decision booleans, decision status, raw process metrics,
artifact hashes, terminal records, and structural counters remain exact.

The screen is not admitted until the corrected runner independently verifies
the preserved `result.json` to a distinct output.  The backend implementation,
frozen thresholds, target, run order, raw outputs, and measurements are not
changed.  Confirmation will use a runner built from the correction commit and
will record that runner separately in provenance.
