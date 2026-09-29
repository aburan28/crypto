# Local ARM64 exactness preflight (outside the promotion measurement)

Two single-process, one-repeat calls ran on `Darwin arm64 T6041` with
Rust `1.98.0` and one Rayon thread, from source commit `32199ca0`.
The same release binary had SHA-256
`ac38adc5e7aa68ce9b84ba754d0835c947b9e585e85acf9f43cc7ab8d91d0ed6`.
Source SHA-256 values were
`cc50fa7dee14d6c606883713fe842139a1d25f636c57be99817ca63607c2b816`
for `koblitz_groebner.rs`,
`e8c66b4d9163160ed32b848820beb5013eb384ddde09a8bc7b4919565ae93959`
for `matrix_f5_f2.rs`, and
`35cafca3dc64e8cf25a1ec63b1338eddb8da4725d0ca1c0f9dedc98313cb16aa`
for `f4_f2_bench.rs`.

Commands (both also set `RAYON_NUM_THREADS=1`, `KIC_F5_ECHELON=2`,
`KIC_F5_FUSED_BUILD=1`, `KIC_GF2_FORCE_AVX2=1`):

```text
KIC_F5_DIRECT_PACK=0 target/release/examples/f4_f2_bench 1 24 f5
KIC_F5_DIRECT_PACK=1 target/release/examples/f4_f2_bench 1 24 f5
```

The raw outputs are [prior-arm64.jsonl](prior-arm64.jsonl) (SHA-256
`f2dfaba5e522ddcf31efa0786951b0e0a2901a7bc1c9cca1891bfe6c1c049372`)
and [direct-arm64.jsonl](direct-arm64.jsonl) (SHA-256
`b5ba3c5d9f1257e57a5022f7ff31df73d36124a1f4246cfff37a4c0b8cf0ccda`).
Every case reported the direct path as used. The prior and direct arms
matched raw row fingerprint, column count and reduction word operations on
all seven cases. This was a preflight for correctness and path activation,
not a paired or calibrated speed measurement. Its one-sample timings are
retained in the raw files but do not enter the frozen Linux gate.
