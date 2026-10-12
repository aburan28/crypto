# Evidence and replay

The first run is `evidence/receipt.json` with SHA-256
`8dfd7c98397db3fa0be5f7b946ae03e831c1942577c758aabb3d67818320fd30`.
It records the exact producer and independent-verifier argv, UTC start/end,
exit, wall, CPU, cumulative child RSS upper reading and SHA-256/byte count
of every child stdout/stderr and JSON artifact. The archived producer result
is `evidence/producer.json` (SHA-256
`dd97bddedae91d41fdedded8d4b6d6b209776952e0c17d4e43771847e0f7609e`);
the independent verifier result is `evidence/verify.json` (SHA-256
`d4519f55a354211a8bdd641d34740e64f7cfb96b8e5c6853a434825ae3183a4c`).

From the repository root, run:

```sh
python3 research/notes/ecc2k130/symbolic_oaware_width_gate_20260925/ci_replay.py
```

That command first checks all frozen source/input bytes, then verifies the
receipt and every committed file hash, recomputes the current width
refusal, and reruns the independent full n13 source/μ4 enumeration from
the archived producer result. The replay's deterministic verifier JSON
must match the committed `verify.json` byte for byte. It does not rerun
the producer or claim an independent n131 solve.
