# Six-summand planted-witness access ladder

This is a frozen SAT-access control for the primary `a=0` Koblitz F₂⁸³ curve. It follows the [capped 4+2 natural-search pair](../hosted-m6-paired-run-39054543374/RESULT.md) on the same sequential K600 base and Gray-prefix or hash-derived K600 residual pools. The three factor-base archives are already stored in the requested S3 prefix. [prepare.py](prepare.py) replays their compressed hashes, full point/phase labels, exact signed-Frobenius overlaps, selected point membership, curve equation, subgroup order, and six-term group sums. The frozen [fixture](fixture.json) has SHA-256 `76b78ec1a68807d130fc11e54e4220b5cde619d677316e313dc247869f727a9d` and BLAKE3 `f07b515f10b4cfd0b0d82267a22c1d46018d278fc922de7cd7b73c77e61ed43a`.

For each arm, the first four summands come from the sequential base and the final two from that arm's exact residual pool. The first pair is ordered by x coordinate. Each arm has its own known-reachable public target, computed from those six subgroup points. It is **not** the frozen random trial used in the natural-search pair, and the two planted targets differ. Thus these runs test accessibility of the same six-summand encoding and domains, not a matched factor-base runtime or natural relation rate.

The ladder uses the same `ChainedS3Encoding`, representative-and-Frobenius-phase domain constraints, first-pair order, zero intermediate-identity mask, 150,000-variable cap, and 12-million-domain-clause cap as the natural search. Its domain-clause counts must equal 237,280 for Gray and 323,796 for hash. For each arm, the stages are:

| Stage | Fixed x blocks | Conflict cap | Outer wall cap |
| --- | ---: | ---: | ---: |
| `all_x` | Six summands and four intermediate sums | 5,000 | 600 s |
| `summand_x` | Six summands; intermediate sums inferred | 10,000 | 900 s |
| `open` | None; only the planted target and full domains remain | 20,000 | 1,800 s |

The two arms run sequentially on one GitHub ARM runner. Each solve is isolated in a no-network, read-only Linux ARM64 container with one CPU, 4 GiB hard memory and swap limits, and 64 PIDs. Before solving, the same pinned Linux image runs the release library and touched-example tests under a 4-CPU, 8-GiB build cap. The job timeout is 180 minutes; conservatively charge the full **10,800-second envelope** including checkout, tests, archive regeneration, build, every attempted stage, failed stages, and artifact upload. The workflow retains the generated archives, manifest, fixture replay, test logs, binary, source hashes, Docker inspection/cleanup, phase markers and raw stdout/stderr. An independent post-run audit must replay those bytes before interpretation.

`SAT_verified_planted_x` demonstrates the solver reached the pinned summand x tuple in a fixed stage after native field-model verification. `SAT_verified_field_model` in the open stage is only a native field-model result; group lifting would be required to call it a partial relation. `UNKNOWN_conflict_cap` and outer wall caps are censored. `PRODUCER_FAILURE_*` marks an encoding, identity, receipt, or infrastructure failure. The planted tuple mathematically satisfies every stage, so an UNSAT solver verdict is a producer failure, not a refutation of a random target. If either arm fails before a stage, later stages are absent and separately classified.

The sweep winner and total index-calculus runtime remain null. Even a positive open-stage planted target does not measure random-target relation probability, large-prime graph closure, useful rank, linear algebra, or an individual logarithm. Its purpose is to distinguish an inaccessible six-summand encoding from an accessible encoding whose natural-search call reached its cap.

To check the fixture against fresh S3 downloads named by compressed BLAKE3 digest:

```sh
PYTHONPATH=PATH_WITH_BLAKE3 python3 research/koblitz_n83_factor_base_sweep_20261008/verification/hosted-m6-witness-ladder-20261017/prepare.py \
  --bases DOWNLOADED_BASE_DIRECTORY --check
```
