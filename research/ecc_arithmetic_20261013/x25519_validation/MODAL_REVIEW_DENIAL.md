# Modal X25519 validation launch was rejected before export

On 2026-10-13, the command below was submitted for automatic approval review:

```sh
modal run research/ecc_arithmetic_20261013/x25519_validation/modal_validation.py --label run1
```

The exact bundle at the rejected attempt is
[MODAL_BUNDLE_rejected.sha256](MODAL_BUNDLE_rejected.sha256):
the new `f25519.rs` and `x25519.rs` sources, five local validation driver
files, and the launcher. The worker request was two vCPUs, 2,048 MiB and a
900-second timeout. The launcher would have run the same correctness checks
on x86-64 and retained its raw outputs.
The historical field source and field driver are preserved as
[`f25519_before_barrier.rs`](f25519_before_barrier.rs) and
[`ecc_f25519_check_rejected.rs`](ecc_f25519_check_rejected.rs).
[`verify_rejected_bundle.py`](verify_rejected_bundle.py) checks all eight
historical hashes without changing the original manifest.

Automatic approval review rejected the action before process creation. Its
stated reason was that general authorization to use Modal for ECC work and
authorization for the earlier seven-file binary-field bundle did not
specifically authorize sending this newly added X25519 source payload to
Modal. It said not to route the payload through an indirect workaround.
No Modal worker or remote X25519 result exists for this attempt.

After the denial, a local ARM64 assembly audit found and fixed two
input-dependent branches in canonical encoding. The current source has a
separately hashed [eight-file bundle](MODAL_BUNDLE.sha256), which has not
been sent remotely. Export remains pending specific user approval of that
current manifest. This is an access condition; it says nothing about
arithmetic correctness or x86 performance.
