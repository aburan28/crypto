# Comparator correction before fresh holdouts

The protocol for this continuation named the five arms of the 2026-09-22
support-envelope study as its reference. That is a reproducible historical
comparison, but it is **not the strongest known constructor for the same
finite matrix workload**. The already merged
`research/boolean_packed_construction_20260922` study retains a packed-direct
arm on n=6/8/10/12, coefficient-changing and degree-cycle families, and
batches including 64. It uses the same ordered Boolean matrix oracle. Its
primary and disjoint confirmation runs report paired 95% lower bounds above
2x in all eight declared batch-64 construction cells against the pointwise
fastest of those five legacy arms. The old source SHA-256 is
`db04d1ae9d9b7fdb8939cd822d48a9a82d637a963d6e3e9d7e1a32631f26061f`;
the primary and confirmation manifest SHA-256 digests are
`af7e15112af3b5136c3a5b281cb4807ebdded2f044f55af2b95e6f363341edaf`
and `8665974a5d8476fcc79cbb95695bff3e9a40e12e3af5398ca0b0e999d6213126`.

Those prior timings used a different host and predate the current isolation
protocol. They cannot be divided into this run's Linux ARM64 timings or
substituted for a same-binary control. They do establish that a future gain
claim must include packed direct in a fresh matched experiment. The frozen
`protocol.json` and both retained attempts here are left unchanged; this is
an additive correction to the scientific interpretation, not a rewritten
gate or a retrospective source change.

This continuation's native discovery is complete and resource-qualified:
128 cells, 8,960 A/B and 2,560 A/A observations, with every output verified.
The indexed arm passed 16/32 incremental groups and **0/32** 2x groups against
the five older controls. The degree-only arm passed 0/32 of either kind. At
n=12 on the two coefficient-changing batch-64 discovery seeds, the pooled
median indexed cost is 5.054149 ms, against 5.789168 ms for the fastest pooled
old control. That descriptive 1.145x ratio is not the paired gate and has no
cross-host comparison with packed direct. The full/holdout gate is unknown.

**Decision:** do not promote this candidate or consume the two fresh holdout
seeds for its declared 2x claim. It missed every dramatic discovery group even
against the weaker historical roster, and its declared reference omitted the
known packed constructor. Future work should start with a matched native
packed-direct control, then test a lever that can beat it. Exact matrix
construction is still only a stage diagnostic; full solving, relation yield,
operation calibration and rho comparison remain unmeasured.
