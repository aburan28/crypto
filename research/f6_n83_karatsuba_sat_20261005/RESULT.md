# n83 F6 recursive carryless circuit: gate target passed, solver timed out

The [preregistered gate](PROTOCOL.md) completed on the same exact
four-link five-summand system as [#1453](https://github.com/aburan28/crypto/pull/1453).
The registered curve was `icv1-f2m83-tm6151469093347-debefd74`; the
standard dimension-18 factor base had 261,447 geometric points, 261,444
subgroup-usable points, and 130,722 sign-folded columns. The ordinary
public T001 and its four rational 4-torsion preimages were unchanged.
The 339 input bits and 332 S3 output equations were unchanged. Each of
the six panel emissions and the additional control emission matched
expanded `System512` on the planted assignment and eight deterministic
arbitrary assignments.

Only the two dense 83×83 products in the middle S3 links changed. The
candidate padded their polynomial operands to 128 bits, recursively
formed three half-size carryless products, reconstructed the 255-bit
unreduced polynomial, and reduced its coefficients under the same field
modulus. Source-coordinate multiplication stayed in the basis-aware
representation. On ordinary offset zero the exact circuit counts were:

| Encoding | AND gates | XOR gates | Mixed constraints | XCNF bytes |
| --- | ---: | ---: | ---: | ---: |
| #1453 basis-aware reference | 30,536 | 68,561 | 160,501 | 2,859,048 |
| Recursive dense multiplier | 20,450 | 50,089 | 111,771 | 2,010,237 |

The candidate had **33.0% fewer AND gates**, **26.9% fewer XOR gates**,
and **30.4% fewer mixed constraints** on this matched input. The other
ordinary offsets had 111,609–112,653 candidate constraints, compared
with 160,339–161,383 in the reference, and all retained 20,450 AND
gates. The preregistered **secondary encoding gate passed**. These are
deterministic representation counts, not measured SAT, F6 or IC speedups.

| Arm | Solver exit | Verification | Peak sampled RSS (KiB) |
| --- | ---: | --- | ---: |
| Planted, all 339 inputs fixed | 10 | `verified_planted`, full-group indices `[0,2,4,6,8]` | 24,768 |
| Planted, only 90 source bits fixed | 124, 30 s | no model | 51,776 |
| Ordinary T001, offset 0 | 124, 120 s | no model | 72,928 |
| Ordinary T001, offset 1 | 124, 120 s | no model | 74,080 |
| Ordinary T001, offset 2 | 124, 120 s | no model | 82,288 |
| Ordinary T001, offset 3 | 124, 120 s | no model | 77,024 |

No process hit the 7-GiB RSS kill gate. All ordinary outcomes are
timeouts, never UNSAT or observed relation misses. The planted witness
contains torsion point index zero, so it is a correctness control rather
than a subgroup-usable IC relation. Thus the primary structural gate
failed: no verified ordinary relation was found under the limits.

The frozen release binary SHA-256 was
`1693d6e5ecf4d70ddbf499aeb7583969b17b8fb37aaa0556fe3645e3b9892293`;
the source SHA-256 was
`de7862943a0adfdb9f8db044f4afe87719befd85f8337ceacd4b924587e34aad`.
The CryptoMiniSat 5.14.7 binary SHA-256 was
`a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac`.
The [runner](run.sh) and [seal script](seal.sh) SHA-256 values were
`62c06828de38a9c6f63c410373c388defefe564b2080874e9385daa8abb828a8`
and `4285f3fcbd68db85e68b83be750ff6328df3649dc5d71809344dd098e57668a5`.
The [status](status.tsv), [RSS samples](rss.tsv), per-arm raw outputs,
independent verifier reports, [input hashes](XCNF_HASHES.tsv),
[planted model hash](SAT_MODEL_HASHES.tsv), [seal log](seal.log), and
[SHA-256 manifest](SHA256SUMS) preserve the complete panel. The
losslessly compressed XCNFs and SAT model were checked against their
original bytes before the originals were removed. The manifest passed
`shasum -a 256 -c` after sealing.

There is no complete ordinary F6 point decomposition, F4/F5 same-input
comparison, one-target IC online interval, verified DLP recovery or
same-point rho reference. The complete candidate ID, operation-count
ratio and end-to-end speedup remain unknown. This contended Apple M4 Pro
has no accepted CPU-isolation receipt; process limits are feasibility
diagnostics only.
