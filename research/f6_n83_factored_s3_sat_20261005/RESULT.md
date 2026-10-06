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
