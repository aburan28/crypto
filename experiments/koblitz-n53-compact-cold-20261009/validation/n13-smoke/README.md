# Compact-orbit replay correctness control

The independent verifier reconstructed all four rank-increasing relations,
four signed-Frobenius base orbits, the solved representative logarithms, and
the target scalar `7` on the exact n13 curve. The strong rho control recovered
the same public point `[384,1476]` and scalar. The archived raw files and their
SHA-256 values are in `receipt.json`. Changing one modular rank-row entry,
one S3 intermediate root, or the reported target scalar made the verifier
reject the copied record.

The producer was built from local cleanup revision `c53f5e7f1` with the
one-line strong-rho CLI fix applied. The receipt pins both source and binary
hashes. This run exercises the replay code and is not part of the n53 timing
panel. Regenerate the archived replay from the crypto repository root with:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 \
  experiments/koblitz-n53-compact-cold-20261009/replay.py \
  --run-dir experiments/koblitz-n53-compact-cold-20261009/validation/n13-smoke \
  --out /private/tmp/n13-compact-replay-check.json
cmp /private/tmp/n13-compact-replay-check.json \
  experiments/koblitz-n53-compact-cold-20261009/validation/n13-smoke/replay.json
```
