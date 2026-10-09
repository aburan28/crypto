# Standard-curve structural comparison (frozen before execution)

Compare SEC 2's `sect113r1` and `sect113r2` on their published binary
field, with `ghs_screen --genus-bound 4` and `ghs_transport` on each published
generator at relative degree 113. Preserve both full JSON responses. These
commands check point membership and annihilation, then the trace-model
precondition. The source constants are `BinaryCurve::sect113r1/r2` in
`src/binary_ecc/curve.rs`, which records the SEC 2 section for each.

For NIST P-192, record its published prime subgroup order and the derived
generic square-root work scale only; the binary GHS command has the wrong
field characteristic and must not be run against P-192. For the existing
`F_(2^192)` fixture, retain its prior structural and trace receipts but
label its order-two subgroup distinctly from a 192-bit prime subgroup.
Neither standard curve is assigned an unmeasured end-to-end index-calculus
cost. The derived sqrt-order scale is a mathematical baseline in group
operations, not a host benchmark.
