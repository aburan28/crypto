# Bounded structural run: 2026-10-12T19:47Z

This is the descriptive run record for [PROTOCOL.md](../../PROTOCOL.md)
and [MAP_PROTOCOL.md](../../MAP_PROTOCOL.md). It contains exact counts,
not a solver or timing result. The [run PDF](report.pdf) embeds the
[diagram](../../fiber.svg) and [census graph](../../census.svg).

- Base commit: `e3090cbc98c657106fef10020da8ccc19be49a3f`.
- Host: local macOS arm64 interactive workspace, no host isolation claimed.
  Compiler: `rustc 1.93.1 (01f6ddf75 2026-02-11) (Homebrew)`.
- Workload: every prime `p≤199`, `p≡1 mod 3`, every nonsquare `b`, plus
  one named map check at `p=13,b=2`.
- Seed: none; enumeration is deterministic.
- Requested and completed: 1,014/1,014 census cells, 1/1 explicit map
  check. Failures, timeouts and invalid rows: zero.

Commands, each exit status 0:

```sh
rustc --edition 2021 -O -D warnings research/isogeny_torsion_fibers_20261012/census.rs -o /private/tmp/isogeny-torsion-fibers-census
/private/tmp/isogeny-torsion-fibers-census > research/isogeny_torsion_fibers_20261012/census.tsv 2> research/isogeny_torsion_fibers_20261012/summary.tsv
cd research/isogeny_torsion_fibers_20261012
rustc --edition 2021 -O -D warnings map_check.rs -o /private/tmp/isogeny-torsion-fibers-map-check
/private/tmp/isogeny-torsion-fibers-map-check > map_check.txt
rustc --edition 2021 -O -D warnings plot.rs -o /private/tmp/isogeny-torsion-fibers-plot
/private/tmp/isogeny-torsion-fibers-plot
cd lean && lake build
```

The raw [census](../../census.tsv) has 349 prime-order rows and the
[summary](../../summary.tsv) gives 1,014 scanned curves across 21 base
primes. Every direct and character count agreed and passed the Hasse
check. `p=13,b=2` has 19 rational points; its degree-3 quotient also has
19. The [map output](../../map_check.txt) reports 19 distinct images and
361/361 ordered-pair homomorphism checks.

Source/output SHA-256:

| file | SHA-256 |
|:--|:--|
| `census.rs` | `226900d93845c885e620b0e4f9f5cfb1efbc6957a7ce7e51c693c24f980927dc` |
| `census.tsv` | `6e36e910172b28dcd8678886b9e94d47785a0c584ab02321a50f48e9138d5dac` |
| `summary.tsv` | `921f2fc1d01a86543bb772f3fa34fcfc9d45dcc2fa32f54cf92ed2fff30b90d5` |
| `map_check.rs` | `6464afc30b059b2016e4e65d6904c998cc6b2b8a6f4e05171c4c0c590156e5a3` |
| `map_check.txt` | `2962c63b518333a3e1efd0f84c2462d412d64a775a4629cd9d787b0bb7b5de3b` |

No cost or uncertainty interval was measured. The graph shows exact
counts by base prime; the source [plot.rs](../../plot.rs) reads only
the recorded summary. There is no ratio against rho in this run.
