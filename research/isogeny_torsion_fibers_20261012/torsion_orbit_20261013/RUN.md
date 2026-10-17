# Native execution receipt

- Protocol commit: `1d01a39a25` (before program execution).
- Host: macOS Darwin arm64, shared/L0. The result is exact support counts;
  wall timing and hardware isolation are not used for a performance claim.
- Runner source SHA-256:
  `8f45a40dc405483a2267bc96f0b48a8f8d839125b0dc8fa8b52420dfb00e54cb`.
- Release native binary SHA-256:
  `2c1193c4ab19b6af615ca6246824d494aa809f3109d29a263a3ff77a6eba3cf8`.
- Raw TSV SHA-256:
  `e4020815166177e7319ecb603f5fb1c32b159307fa95a6bfd729ca19c3f6926d`.
- Retained first compile failure: [compile-v1.log](compile-v1.log), with its
  pre-correction source hash. No experiment ran before the fix.

Commands, from repository root:

```sh
rustc --edition=2021 -O research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/screen.rs -o /private/tmp/isogeny-torsion-orbit-screen-20261013
/private/tmp/isogeny-torsion-orbit-screen-20261013 > research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/raw-v1.tsv 2> research/isogeny_torsion_fibers_20261012/torsion_orbit_20261013/run-v1.log
```

Both commands exited zero after the retained first compile failure. The
output has 12 rows: three arms on each of four curves. All asserted
homomorphism, torsion, curve-equation, scalar-index, subgroup-order and
point-sum/scalar-sum checks passed. The [run log](run-v1.log) has one
verification line per curve. `screen.rs` contains two arithmetic views,
point addition and a full subgroup scalar index, but these are in one
executable. A second same-host execution produced
[replay-v1.tsv](replay-v1.tsv) and [replay-v1.log](replay-v1.log); `cmp`
returned zero for both raw and log pairs. Independent-host replay has
not been performed.
