# Validation record

The first supervisor test compile failed because the included repository SHA256
tests require the `hex` crate. Added `hex` as a development dependency; this does
not change the launch executable's dependencies or map acceptance.

Original diagnostic:

```text
error[E0433]: failed to resolve: use of unresolved module or unlinked crate `hex`
  --> ../../src/hash/sha256.rs:342:9
342 | hex::decode(hex_str).unwrap()
error: could not compile `additional-prime-isogeny-campaign` (bin "isogeny-campaign" test)
exit status: 101
```

The updated CLI passes 63 tests plus 3 field-reference tests. The native
supervisor passes 14 tests, including normalized map-identity checks and the
repository SHA256 reference vectors. The root release library gate and walk
replay controls are running for this source revision. The prior merged-tree
validation remains historical evidence and does not replace the new checks.

The first commit included the CLI changes but sparse-checkout staging omitted
the new study directory. The following protocol commit explicitly stages the
complete study and lockfile before any measured campaign execution. The clean
producer revision is recorded in LAUNCH.json.
