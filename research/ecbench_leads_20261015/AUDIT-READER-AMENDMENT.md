# Extend the v1 decimal reader without changing recorded evidence

The four exact E0 N83 attempts executed and exhausted their frozen query caps.
The first audit/table/resource commands then failed before reading the records:
`1.5030185107919957e-7 is outside the v1 decimal-parser domain`.
The command manifest and original failure logs are preserved. No sealed session,
method, budget, count, scalar, derived figure or record seal is edited.

The compatibility helper had only23 powers of ten, through10^22. The historical
serde_json non-float_roundtrip parser it reproduces has309 literal powers,
through10^308, and scales very small exponents in chunks. A17-digit significand
for this tiny S needs10^23. This is an implementation-domain gap, not a numerical
tolerance failure or a cryptanalytic obstruction.

Extend the helper to the historical parser's complete power table and underflow
scaling. Preserve every original table value and arithmetic path for exponent
magnitudes0..22. Do not weaken integer comparisons, seals or replay checks.
Validate against an independently pinned historical serde_json reader, existing
bit-pattern tests, finite extremes and the exact preserved N83 session. Store a
new audit receipt and separate resource/table outputs with the corrected reader;
replay unsuccessful executions with the original measured binary. Record auditor
and measured binary identities separately. Original successful audits stay intact.

The registered IC algorithm is unchanged. A tiny counted S on an exhausted N83
attempt is censored work, not an attack cost or a win against rho.

The first full-library check retained4094 passes and one new reference-test
failure. That test incorrectly identified raw-token parsing with the existing
serialize/parse contract: `1.2345678901234567e-9` has a different canonical
serialization, and the historical integer-significand parser can distinguish
the two by one bit. Factor the raw-token parser from the serialize/parse wrapper.
Check all14 unchanged historical reference bits on raw tokens and additionally
check the exact actual N83 token through the replay roundtrip. No expected bit
pattern or tolerance is changed. Preserve the failed library output separately
and rerun the required checks.
