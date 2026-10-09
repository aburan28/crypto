# Compact-orbit replay correctness control

The native independent verifier reconstructed all four rank-increasing relations,
four signed-Frobenius base orbits, the solved representative logarithms, and
the target scalar `7` on the exact n13 curve. The strong rho control recovered
the same public point `[384,1476]` and scalar. The archived raw files and their
SHA-256 values are in `receipt.json`. Its replay JSON matches the archived
output byte for byte. Rust tests change one modular rank-row entry, one S3
intermediate root, or the reported target scalar and require rejection.

The producer was built from local cleanup revision `c53f5e7f1` with the
one-line strong-rho CLI fix applied. The receipt pins both source and binary
hashes. Its original verifier used Python; those historical source and result
files remain archived. This run exercises the native replay code and is not
part of the n53 timing panel. Regenerate its replay from the crypto repository
root with:

```sh
CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/rho-strong-cli-20261009/target \
  cargo build --release --example koblitz_n53_compact_replay
/Volumes/SSD990/crypto/worktrees/rho-strong-cli-20261009/target/release/examples/koblitz_n53_compact_replay \
  --run-dir experiments/koblitz-n53-compact-cold-20261009/validation/n13-smoke \
  --out /private/tmp/n13-compact-replay-check.json
cmp /private/tmp/n13-compact-replay-check.json \
  experiments/koblitz-n53-compact-cold-20261009/validation/n13-smoke/replay.json
```

The `CARGO_TARGET_DIR` above reuses this machine's build cache; another host
may use any writable Cargo target directory and its corresponding binary path.
