# Required check repairs

The first PR head `54eff85de762435da4456bbaf0bb9ce846705fe5` failed Rust
formatting and one historical-source test. These failures are retained in the
[lint run](https://github.com/aburan28/crypto/actions/runs/36823472116),
[autolab run](https://github.com/aburan28/crypto/actions/runs/36823472003) and
[producer run](https://github.com/aburan28/crypto/actions/runs/36823472001).
No registered experiment was dispatched by these pull-request checks.

The serializer's longer JSON expression requires an `Ok(...)` line wrap under
the pinned formatter. Applying that formatting changes source bytes, not the
report contract or algorithm. The local native-control archive retains its
original pre-format source binding; CI rebuilds the final formatted source.

`test_mathematical_fixture_is_shared_and_contains_no_preparation_history`
incorrectly compared today's live worker to the immutable v1 source pin.
The repaired test hashes the actual worker in the retained v1 source archive.
It continues checking the shared mathematical fixture and lack of historical
preparation inputs. A new independent negative test constructs a changed source
archive with internally consistent manifest, build policy and identity digests;
the production v1 admission gate must still reject it specifically at its
unchanged worker-source pin. No asset, production validator, source pin, old
registration or raw experiment was rewritten or weakened.

The targeted module's ten tests pass locally, including the retained-source
check and consistent-metadata rejection. All applicable checks must pass again
on the final PR head; earlier successes do not substitute for that validation.
