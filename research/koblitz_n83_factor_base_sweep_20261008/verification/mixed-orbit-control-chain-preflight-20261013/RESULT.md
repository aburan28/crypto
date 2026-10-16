# Known-answer N83 chained-S3 preflight

The frozen five-point K=1,182/K=2,048 mixed-orbit control fixture has a
concrete witness for the **all-affine identity pattern** of the worker's
five-summand chained S3 sieve. [The replay](replay.json) records zero for all
four S3 polynomial residues and records the three intermediate x-coordinates
and target x-coordinate. The first two small-domain points must be swapped to
meet the worker's `x0 <= x1` ordering; their interchange preserves the group
sum. The ordered global point indices are `166, 0, 332, 196212, 196378`.

[The checker](replay.py) pins the exact fixture SHA-256
`fc203db28eed0164e47b0c755c3d9a09f80df1acd5e4c7e09dfe172fa5f0a538`
and its two retained compressed-object hashes. It uses a separate carryless
polynomial multiplier and generic binary-curve addition to recompute the
target and each S3 value over `F2^83`. It checks that each selected point is
on the curve, that all terms are the phase-zero representatives with the
declared small/residual domains, and that the full ordered group sum matches
the fixture target. The checker also pins the audited K=600 gate and byte
hashes for the exact three worker/chain source files from frozen build commit
`f3d27f68b68cda907e0668042ac7fa13ae6b579d`. The earlier
[archive-bound fixture construction](../mixed-orbit-known-witness-source-20261012/RESULT.md)
is the evidence that those selected terms came from the exact retained S3
objects. This replay does not reopen the compressed archives.

Three [mutation checks](test_replay.py) pass: a changed fixture byte is
rejected, an incorrect S3 output has a nonzero residue, and reversing the
canonical first pair is rejected. The saved receipt reproduces byte for byte:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-control-chain-preflight-20261013/replay.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/mixed-orbit-control-chain-preflight-20261013 -p 'test_*.py' -v
```

This proves the field-level necessary conditions and the selected group
equation for one constructed witness. It does not build a SAT assignment,
measure whether CDCL reaches that witness, or establish natural relation
yield. The next [bounded solver-access control](../mixed-orbit-control-reuse-17400-20261012/RESULT.md)
still needs its separate 1,200-second authorization. This source-only check
launched no Docker container or SAT solver and made no new local-budget charge;
the conservative ledger remains 17,400 of 17,400 authorized seconds.
