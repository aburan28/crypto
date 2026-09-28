# ECC2K-130 degree-263 isogeny preflight

**Observation, exact public target.** For ECC2K-130, the saved deterministic
preflight constructs 264 degree-263 kernel lines from rational torsion on the
quadratic twist. The quotient parameters give two horizontal loops and 262
distinct descending j-invariants, arranged in two Frobenius orbits of 131.
The independent replay checks oriented full-point images for one kernel from
each of four line orbits in each of two seeded runs. This establishes a
reproducible representation construction and representative map checks. It
does not establish cheaper point decomposition, useful relation yield, or an
ECDLP speed improvement.

The pre-registered boundary is [twist_torsion_contract.md](twist_torsion_contract.md).
[twist_torsion_preflight.py](twist_torsion_preflight.py) imports the tracked
`research/ecc2k130_relations/fastfield.py` and `relations.py`, then records
their SHA-256 hashes. [twist_torsion_results.json](twist_torsion_results.json)
contains both seeds and the representative kernels. The independent
[redteam_velu_replay.py](redteam_velu_replay.py) imports no production
finite-field or group arithmetic; its [result](redteam_velu_replay.json)
reports 8 representative kernels, 32 codomain-valid point images, 16 full
additivity equalities, and zero disagreements. The formula and remaining
limits are explained in [isogeny_theory.md](isogeny_theory.md).

From the repository root, reproduce without replacing the archived evidence:

```sh
python3 research/ecc2k130_direction_review_20260924/twist_torsion_preflight.py --out /tmp/ecc2k130_twist_preflight_replay.json
python3 research/ecc2k130_direction_review_20260924/redteam_velu_replay.py --input /tmp/ecc2k130_twist_preflight_replay.json --out /tmp/ecc2k130_velu_replay.json
```

The submitted preflight was freshly rerun on both frozen seeds (3.27 and
3.78 seconds); the independent replay then passed all 8 representative
kernels in 8.40 seconds. The deterministic payloads and diagnostic counters
of both new JSON files agree exactly with the earlier saved runs after
excluding timestamps, elapsed time, command paths and hash metadata. The
source, contract, dependency and input hashes in each JSON artifact identify
the files used in the submitted run. Timings and operation counters are
diagnostic, not calibrated comparisons. The replay covers only saved nonexceptional points;
infinity and kernel inputs, a generic production y-map, dual composition,
factor-base transport, relation collection, and linear algebra remain open.
