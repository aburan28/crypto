# Preserve search completeness while testing faster coordinates

Registered before changing the search status code or running performance
measurements. The first candidate release suite failed for two reasons:
the new test used an unavailable `KoblitzCurve::new(0,17)` fixture, and
the existing suffix-recursion test intermittently stopped when a
parallel branch encountered an `Infinity` remainder. The fixture is
corrected to the available `KoblitzCurve::new(1,17)`. The initial full
test failure is preserved as `wide_tests_initial_failure.log.gz`.
After the fixture fix, exact coordinate lists matched the old curve-lift
reference for n=9, 17 and 83, restricted indices and a nonstandard
basis. The suffix test still failed in one of ten direct replays.

The suffix bug predates indexed coordinate screening: internally,
`unsupported` conditions also set the global `exhausted` flag, causing
one unsupported branch to stop other branches that might contain a
witness. Change internal unsupported branches to set only
`stats.unsupported`; continue searching other branches. At the public
`wide_groebner_decompose` boundary, a no-witness result with any
unsupported branch must still set `exhausted=true`, preserving callers
that use that flag to avoid a false proof of absence. A verified
witness clears the incomplete flags. This is an exactness correction,
not a performance hypothesis. The n83 root probe should never enter an
unsupported branch, so the paired root timing remains a controlled
comparison of coordinate screening.

Before performance measurement, require all wide tests to pass and
20 direct suffix-test replays to pass with the candidate binary.
Preserve all failures and report any continued flakiness; if the
correctness gate fails, do not promote the coordinate optimization.
