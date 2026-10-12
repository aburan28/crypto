# n83 F6 basis-aware S3 SAT encoding: smaller input, no ordinary model

The [preregistered gate](PROTOCOL.md) completed on the exact four-link,
five-summand system inherited from [#1452](https://github.com/aburan28/crypto/pull/1452).
The curve was `icv1-f2m83-tm6151469093347-debefd74`, with a dimension-18
standard base of 261,447 geometric points, 261,444 subgroup-usable points
and 130,722 sign-folded columns. The same public T001 target and all four
rational 4-torsion preimages were used. The emitter kept the same 339
inputs and 332 S3 output equations. Every one of the six panel inputs and
the separate control emission matched the expanded polynomial system on
the planted assignment and eight deterministic arbitrary assignments.

The characteristic-two identity
`S3(x,y,z) = (xy + xz + yz)^2 + xyz + b` was implemented with the 18-bit
source coordinates represented directly in their field basis. On ordinary
offset zero, compared with the prior native-XOR circuit, AND gates fell
from 36,761 to 30,536 (**16.9%**), XOR gates from 106,013 to 68,561
(**35.3%**), and mixed constraints from 216,628 to 160,501 (**25.9%**).
The other ordinary offsets had 160,339–161,383 candidate constraints
versus 215,304–217,474 reference constraints. The reference and candidate
are different circuit representations of the same equations, so these are
exact encoding counts, **not** measured solver or F6 speedups. The
preregistered secondary gate required a 20% drop in both AND gates and
constraints; it failed on AND gates.

| Arm | SAT exit | Independent verification | Maximum sampled RSS |
| --- | ---: | --- | ---: |
| Planted, all 339 inputs fixed | 10 | `verified_planted`, full-group indices `[0,2,4,6,8]` | 28,848 KiB |
| Planted, only 90 source bits fixed | 124, 30 s | no model | 62,608 KiB |
| Ordinary T001, offset 0 | 124, 120 s | no model | 84,848 KiB |
| Ordinary T001, offset 1 | 124, 120 s | no model | 85,536 KiB |
| Ordinary T001, offset 2 | 124, 120 s | no model | 83,328 KiB |
| Ordinary T001, offset 3 | 124, 120 s | no model | 80,384 KiB |

The single maximum across the runner's sampled RSS rows was 85,536 KiB;
no process triggered the 7-GiB kill gate. All four ordinary outcomes are
timeouts, not UNSAT or observed relation misses. The planted witness uses
torsion point index zero and has no subgroup-usable witness, as in #1452;
it is a correctness control only. The [status](status.tsv),
[RSS samples](rss.tsv), per-arm raw outputs and verifier reports preserve
every attempted case. [XCNF_HASHES.tsv](XCNF_HASHES.tsv) gives both
uncompressed and lossless gzip hashes for every generated input, including
the additional pre-panel control. [SAT_MODEL_HASHES.tsv](SAT_MODEL_HASHES.tsv)
does the same for the full planted SAT model. The first [seal log](seal.log)
records a failed checksum caused by hashing that log while it was still
being written; a fresh post-exit [SHA-256 manifest](SHA256SUMS) was verified
successfully. No arithmetic or solver result changed.

The frozen release binary SHA-256 was
`ccf7fdedf814b2464b9f920079d04d8bfb1648be3feddc18fb8478793d122e01`;
the source SHA-256 was
`ada17b8a1d4fb04f06b7038e90d156146f3460867c45c9f455b08681a91d2f74`.
The CryptoMiniSat 5.14.7 binary SHA-256 remained
`a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac`.
The [runner](run.sh) SHA-256 was
`9a750ecb652f3e22ece0dabb291e5a8cb04ed7970c7c65bae9e24238dbd3d99c`.

No ordinary relation, complete F6 point decomposition, one-target IC
online interval, recovered scalar or same-point rho reference was measured.
The candidate ID, operation-count ratio and end-to-end speedup remain
unknown. The contended Apple M4 Pro host lacks an accepted CPU-isolation
receipt, so these process limits are feasibility diagnostics only.
