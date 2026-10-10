# n83 factored S3 circuit: CNF is compact, planted Kissat solve times out

The [preregistered initial gate](PROTOCOL.md) constructed the exact
five-summand S3 circuit for K0 curve
`icv1-f2m83-tm6151469093347-debefd74`. Each of its 332 output bits
matched the existing expanded `System512` polynomial on the planted
full-group assignment and eight deterministic assignments. This is an
algebraic encoding check; it does not show that SAT can recover an
unfixed point decomposition.

The pinned planted target was the full sum of source factor-base indices
`[0,2,4,6,8]`. Factoring the four S3 links gave 36,761 AND gates and
105,060 XOR gates; standard Tseitin CNF had 142,160 variables,
530,945 clauses, 1,518,469 literals and 422 unit clauses, including
the 90 fixed source bits. The CNF was 11,148,377 bytes, SHA-256
`ce469726b1800f6a0cf113d341429a0effbf8019db3e416876b21dde384bc59d`.
Its [lossless gzip copy](planted_0.cnf.gz) has SHA-256
`8e1db554f464888b785f0db59e37bd9fbef4be6a6908271f23fc52a3e44d0dc7`;
decompression compared byte-for-byte with the original. The expanded
Boolean system had 1,091,334 term occurrences. These are different
representation units and do not constitute a speed ratio.

Kissat 4.0.4, binary SHA-256
`05d6f3e9c402a1fe8853b0746e384e1b3d1c4a550e255f11daa2461d279aa848`,
exited 124 at the registered 30-second planted limit without a model.
Its own log records 157,941 conflicts and 858,016 decisions before
termination. The live `ps` sampler saw at most 72,960 KiB resident
before the solver changed its allocation pattern, and no 7-GiB kill.
Kissat's printed macOS resident-size conversion is inconsistent with
the host samples, so it is not used as the measured peak. The verifier
reported `no_model`; the runner correctly stopped before ordinary T001
offsets. This is an inconclusive solve, not UNSAT and not a verified
decomposition.

The native probe source SHA-256 was
`7db5030cc8322b59bc1618a04db5754267f87f84b2738caecd94f79487c88a81`;
release binary SHA-256 was
`a8ef4cf4c9f6fd87026295df7aed87bac12c33fec1580cbc9eb67c969f055b94`.
The first build's borrow error and successful corrected build are both
retained. Raw [solver output](planted_0.solver.stdout), stderr,
[status](status.tsv), [RSS samples](rss.tsv), emission receipt, and
[runner](run.sh) preserve the attempt. No natural relation, complete
F6 call, one-target IC interval, paired rho reference or speedup was
measured. An amended XOR-native solver gate must be registered before
new measurements.

## Amendment 1: native XOR rows

[Amendment 1](AMENDMENT_1.md) was committed before changing the emitter.
The exact same field circuit was encoded with one native XOR row per XOR
gate and three CNF clauses per AND gate. For the planted target, this
reduced the stated constraint count from 530,945 CNF clauses to
216,014 mixed constraints. The two counts represent different solver
languages; their ratio is a representation diagnostic, not a runtime
speedup. Each of the six amended inputs passed the nine-case comparison
against the expanded 332-equation system. [XCNF_HASHES.tsv](XCNF_HASHES.tsv)
records the uncompressed byte counts and hashes of all losslessly
compressed inputs.

CryptoMiniSat 5.14.7, binary SHA-256
`a3f85c3709b5e2a040bf82a4a604d1c7b9f10219bbf180a9e0f72319a2e892ac`,
returned SAT for the fully pinned planted assignment. The model made
all expanded equations zero and replayed the exact full-group source
indices `[0,2,4,6,8]`. Point index 0 is torsion and its cofactor-four
image is the identity, so this is **only a correctness control**, not a
usable IC relation. The first verifier erroneously classified it as an
algebraic-only model by requiring subgroup usability in the planted
control. Its raw model and report remain intact; the corrected verifier
replayed that same model and reported `verified_planted` before the
panel was rerun. Ordinary relations still require every source point
to have nonidentity subgroup image.

| Arm | Solver exit | Verification | Maximum sampled RSS (KiB) |
| --- | ---: | --- | ---: |
| Planted, all 339 inputs fixed | 10 (SAT) | full-group replay; unusable for IC | 30,656 |
| Planted, only 90 source bits fixed | 124 (30 s limit) | no model | 77,856 |
| Ordinary T001, torsion offset 0 | 124 (120 s limit) | no model | 98,432 |
| Ordinary T001, torsion offset 1 | 124 (120 s limit) | no model | 94,832 |
| Ordinary T001, torsion offset 2 | 124 (120 s limit) | no model | 97,344 |
| Ordinary T001, torsion offset 3 | 124 (120 s limit) | no model | 95,056 |

All six processes stayed below the 7-GiB observed RSS gate; none was
killed for memory. The ordinary extended-DIMACS files each had
36,761 AND gates, 104,689–106,859 XOR gates and 215,304–217,474
mixed constraints. Every ordinary run timed out, so no absence of a
five-summand solution is established. The corrected probe source
SHA-256 was
`ae34ec65ff718bff41819fc0f6b206bf71309a1d83180af9675d03e2dd45183e`;
its frozen release binary SHA-256 was
`2261bd409f6754cb23242a2dfcda20e8a1f2a5d3efbb2873e2032fa6b848f619`.
The [XOR runner](run_xor.sh), [status](xor_status.tsv), [RSS samples](xor_rss.tsv),
raw solver outputs, verifier receipts and [SHA-256 manifest](SHA256SUMS)
preserve the panel. Native XOR representation alone did not complete
an ordinary n83 decomposition under these limits; F6 and end-to-end IC
speedups remain unknown.

The two 1,041,420-byte planted SAT model outputs are retained as
byte-verified gzip files; [SAT_MODEL_HASHES.tsv](SAT_MODEL_HASHES.tsv)
records each original SHA-256 and compressed SHA-256. The ordinary
timeout outputs remain as their exact original files.
