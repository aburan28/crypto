# Independent replay of the retained n=53 compact-orbit relations

`replay.py` reruns the twelve public synthetic fixture labels in the archived
compact-orbit panel (four planted, eight natural). It pins the original Rust
source and executable hashes, the certified point-set hash, the original
results hash, and the selected JSONL base header. Each subprocess uses the
orbit-factorized configuration from the IC boundary skill and a fixed 2,000
conflict limit. The Rust SAT solver checks the extracted, pinned witness; it
is not performing an unrestricted relation search in this arm.

The separate Python verifier implements carryless polynomial convolution and
long reduction, exponentiation-based inversion, direct affine group addition,
and the S3 polynomial. It checks all 23,320 retained base points are on the
curve, recomputes natural targets as published scalar multiples of the
generator, and recomputes planted targets from their published point indices.
For each archived and newly extracted x tuple, it checks that a sign choice
uses retained subgroup points, matches **both** pinned pair-sum abscissae, and
sums to the public target. The three S3 identities of the balanced chain are
also recomputed in the independent field arithmetic.

The JSONL container available for replay has the same certified base hash
`d859319…` as the archived run, but a different whole-file BLAKE3 because
the retained files append different observations after the base header.
Both fingerprints are recorded in `validation.json`; this is a point-set
replay, not a byte-identical input-file replay. The same pinned executable
regenerates fixture scalars and planted indices. HashMap iteration can pick
different valid relations on replay; the archived relations are verified
separately rather than requiring identical witness choices.

`run_003/validation.json` is the result, with the raw producer JSONL and
stderr for all twelve runs. A passing result validates these fixed relations
and the observed 12/12 replay, including 8/8 natural. It does not measure a
full-rank collection, a final logarithm solve, a `vs_rho` crossover, or an
asymptotic speedup. The official boundary ledger was not promoted.

Reproduce from the `/Volumes/SSD990/crypto` repository root using a new
output directory:

```sh
python3 research/sat_factor_base_review_20260908/autolab_orbit_extract_20260924/independent_replay_20260924_codex/replay.py \
  --out /tmp/kic-orbit-replay-new
```
