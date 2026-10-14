# N83 known-witness algebra replay of the first identity pattern

The [replay receipt](replay.json) confirms that the pinned public v7 five-point
known-answer fixture has an affine partial sum after every term. Its final sum
equals the pinned target, and the four binary Semaev S3 residuals are all zero.
The native mixed-orbit worker visits `identity_mask=0` first, so this witness
has the algebraic shape of that first all-affine pattern. The v7 worker spent
its 5,000-conflict cap on one identity pattern and returned
`UNKNOWN_search_cap` with zero SAT models. That bounded result does not refute
the known group/S3 witness.

`replay.py` reads only committed files. It pins the fixture, v7 result and
independent audit, Python reference arithmetic, and both relevant Rust source
files by SHA-256, then checks each against baseline commit
`d42167ff4e80a9507911338df35572b5d24021c5`. The fixture uses three
small-domain and two residual-domain positive phase-zero points at global
indices 0, 166, 332, 196212, and 196378. The replay computes each partial
sum with the existing Python group implementation and evaluates
`S3(x,y,z) = xyz + (xy + xz + yz)^2 + 1` in the same binary field. It checks
that all four equations vanish and the final affine point matches the target.

Reproduce from the repository root with:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/known-witness-first-pattern-replay-20261014/replay.py --check
```

This is a source-only algebra check. It does not build or solve the SAT CNF,
prove the complete encoding satisfiable, sample natural relations, or add any
N83 search time to the fully charged 21,000-second ledger. It therefore does
not select a factor base or establish a complete cold one-target runtime.
The frozen [v7 control](../mixed-orbit-control-reuse-19800-20261013/RESULT.md)
and its independent audit retain the empirical solver outcome.
