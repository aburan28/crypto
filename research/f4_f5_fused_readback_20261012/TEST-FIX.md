# Verification repair found by the mandatory library check

The first corrected-candidate full-library run failed: 3,994 passed,
one failed and 112 were ignored (303.77 seconds, exit 101). Its
[log](records/attempt-2-fixed/library.log), catalogue, binary hash and
[single-test reproducer](records/attempt-2-fixed/coordinate-fixture-reproducer.log)
are retained. This is a real failed validation run, not a passing check.

The failure is
`curve_traits::tests::unavailable_coordinate_models_keep_exact_fields_without_order_certificates`.
The inline Montgomery fixture over field size 101 incorrectly supplied 101
as its curve order and expected `NotEvaluated`; the now-supported conversion
and exhaustive counter reject that order (count 96). The unavailable-extension
fixture still correctly has `NotEvaluated` status. The coordinate module and
checker sources are identical to the pre-optimization baseline, as recorded
by the [empty baseline diff](records/attempt-2-fixed/coordinate-baseline-source-diff.txt).
Neither calls the modified F4/F5 code in this test.

Scope expanded to the test module in `src/cryptanalysis/curve_traits/mod.rs`;
the unavailable [Conductor reservation](records/fixture-conductor.log) is
retained. Correct the inline test, not production curve arithmetic or a
registered curve: keep both exact-field checks, require the supported
Montgomery order to be `Proved` without a generator certificate, explicitly
assert that the old wrong order is rejected, and keep the unavailable-extension
assertion. Rename the test to describe both behaviors. Do not disable or
weaken the checker, change a registry/ledger entry or count the old run as passed.

Rerun the affected coordinate tests and then the complete locked release-library
check with the same three test threads and 20-minute budget. Keep fresh logs,
catalogue and exit status under `records/verification-final/`; leave the failed
run unchanged. Verify Cargo example targets too. The repair is inside
`#[cfg(test)]`; the measured production F4/F5 source and comparison remain
unchanged. Record the rebuilt driver's hash if Cargo rebuilds normal artifacts.

Current main has 112 ignored library tests, including the previously long
non-weak isogeny test (the old hashing study had 110 ignored). Keep those
upstream ignore attributes unchanged and disclose the standard-suite coverage.
Run that exact unmodified isogeny test explicitly with `--ignored` in a separate
copied release-test binary, normal Rayon defaults and its own 20-minute budget.
The copied binary precedes this independent fixture-only repair; its isogeny
test/source are unchanged. Retain its exact command, hash, timestamps and
result separately. A standard-suite ignore is not a pass of that test.

The first fixture-repair rerun passed seven coordinate tests and failed the
renamed test on its previously unreached extension input: `missing coefficient
array modulus`. That failed log remains under `records/verification-final/`.
Complete the synthetic extension fixture with the required canonical
`modulus:["2","0"]` (implicit monic leading coefficient); keep its exact
field-size/NotEvaluated assertions and explicitly reject a record with the
modulus removed. No parser rule is relaxed. Fresh final verification records
use `records/verification-final-2/`. This corrects the earlier statement that
the extension fixture was already valid: its intended unavailable-order
status was unchanged, but its serialization was incomplete.
