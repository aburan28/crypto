# Implementation evidence

Classification: accounting and reproducibility. This changes the launch/audit
contract, not an IC algorithm or its measured speed. No new production
registration, native solver invocation, fresh target or paired measurement was
created. All scientific speedup and promotion fields remain unset or false.

The macOS ARM64 development doctor found the source and native build tools,
but no Valgrind; the calibrated Linux instruction protocol is unavailable on
this local host. The following Python 3.12 correctness controls passed; their
logs and content digests are retained in [local-validation.json](local-validation.json).
These suites overlap and their counts must not be added as independent tests.

| Check | Result | Scope |
| --- | --- | --- |
| Prepared target, SAT, F5 v1/v2, shared transport and launch regression suite | 47 passed | Before the final additional rejection/relocation cases |
| Final shared transport, SAT and F5 v2 suite | 20 passed | Isolated transport, rejection, timeout and existing family controls |
| Final claim suite and three real-family result-contract controls | 15 passed | Claim race, distinct-output retries, relocation and auditor field compatibility |
| Skill metadata validation | Passed | Repository-local IC skill |
| `git diff --check` | Passed | Whitespace/diff integrity |

The first 17-test development run passed all 11 claim controls and failed six
transport fixtures before their synthetic entrypoint could run: the fixture's
empty `producer` directory had no source to retain. Adding `producer/control.py`
corrected the fixture. No package source gate was relaxed. The corrected runs
above retain the failure's cause rather than treating that first run as a pass.

The new transport tests execute the actual public CLI and a fresh isolated
frozen helper. The known-name mathematical auditor in those small fixtures is
synthetic and explicitly has no verified scalar or online cost. Real family
result compatibility is separately exercised with the existing preparation
and witness controls, whose native/source calls are mocked. Neither kind of
test constitutes an actual prepared native execution or measured natural yield.

Applicable PR CI is required on the final head before merge. Measured workflow
jobs require their separate dispatch conditions; PR correctness controls do not
reopen any historical confirmation or consumed registration.
