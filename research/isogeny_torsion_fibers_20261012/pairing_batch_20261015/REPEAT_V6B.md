# V6B scheduled same-input repeat after release tests

Frozen before the repeat. The [V6 spec](RUN_SPEC_V6.json), native source
SHA-256 `ad0cab1eeab279f3f47344b2b94e5f93998861df0e4bde78ebede904f1d37015`,
public order-67 target list, H selection, three arms, prefix list, round
indices 0–2, arm orders, rho cap and retry bound are unchanged. Run the
same compiled release binary after `cargo test --release --lib` has exited.
Keep the V6 run and its concurrent-test context intact. Write V6B to new
`*_V6B` files and independently replay all 594 scalar outputs.

Hypothesis: the direct `F_401` negation-rho arm remains below the pairing
arm at `N=64` and `N=66` when the concurrent release test suite is absent.
Reference: direct `F_401` rho, same targets, same public-input charge.
Success for a comparative observation requires all three arms to solve the
same targets, replay 594/594 scalars, and report complete per-prefix cold
current-thread CPU, wall, attempts, primitive counts and RSS. An exhausted
rho or replay mismatch censors the affected round. Stop after exactly three
rounds and one independent replay; retain every failure. This is still a
shared-host L0 repeat; no isolated wall-time or formal counted-`S` claim.

Run commands, with the same paths and cell as V6:

```sh
target/release/examples/torsion_pairing_batch compact-prime-run \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/INPUTS_V4.jsonl 2 \
  "$round" research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/PUBLIC_401.json \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/SELECTION_401_V1.json
target/release/examples/torsion_pairing_batch compact-verify \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/INPUTS_V4.jsonl 2 \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/PUBLIC_401.json \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/ANSWERS_401.json \
  research/isogeny_torsion_fibers_20261012/pairing_batch_20261015/COMPACT_PRIME_401_V6B.jsonl
```

The initial pre-run draft contained mistaken V4 suffixes; these were
corrected before execution. It also used cell index zero. The first
invocation failed its `p=401` assertion before timing or scalar output;
its stderr is preserved as `COMPACT_PRIME_401_R0_V6B_INDEX0_FAILURE.stderr`.
`INPUTS_V4.jsonl` contains `(61,19)`, `(271,61)`, `(401,67)` in index
order, so this corrected command uses cell two. The failed invocation
does not count as a measured round.
