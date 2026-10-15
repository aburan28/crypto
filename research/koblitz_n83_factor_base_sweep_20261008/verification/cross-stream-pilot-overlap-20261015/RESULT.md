# Exact primary cross-stream signed-orbit overlaps

This is an **independent replication** of the earlier, more complete [2026-10-13 exact overlap census](../pilot-overlap-census-20261013/RESULT.md), which had already measured all **648 ordered distinct-stream pairs** among the 27 primary-arm pilot archives and checked binary-curve arithmetic. The new [receipt](overlap.json) reparses the archived point records from the later S3 replay cache: nine `(policy, seed)` streams, each at 64, 256 and 600 signed-Frobenius orbit columns. Each local compressed archive matched its committed pilot manifest and the prior 57-object S3 content audit by BLAKE3, and its decompressed bytes matched the manifest. The generator checked all 166 point records, 83 phases, two signs and Frobenius coefficients in every orbit before computing point-set-derived orbit intersections. It also retained each compressed object's SHA-256.

The [new independent audit](audit.json) reparsed all 27 archives and reconstructed each complete orbit with a different identifier: its lexicographically least point. It checked representative rows, complete phase/sign labels, disjointness within each archive, the S3/manifest binding and all 648 ordered intersections and residual counts. The [comparison receipt](comparison.json) then matched **every one of the 27 archive identities and 648 ordered overlap/residual rows** to the 2026-10-13 census. Two focused corruption controls reject a duplicate phase and a changed residual count. This extra replay corroborates the existing result; it does not add a new solver measurement.

| Pair class, either orientation | Ordered pairs | Disjoint orbit sets | Empty pool-minus-base residual | Shared-orbit range |
| --- | ---: | ---: | ---: | ---: |
| Gray prefix with Gray prefix, different seeds | 54 | 0 | 6 | 63–599 |
| Gray prefix with hash | 162 | 162 | 0 | 0 |
| Gray prefix with sequential | 162 | 0 | 18 | 40–504 |
| Hash with hash, different seeds | 54 | 54 | 0 | 0 |
| Hash with sequential | 162 | 162 | 0 | 0 |
| Sequential with sequential, different seeds | 54 | 0 | 9 | 62–599 |
| **All ordered pairs** | **648** | **378** | **33** | **0–599** |

Thus **615** pairs have a nonempty exact residual signed-orbit domain. Of these, **57** have exactly one residual orbit; the current mixed-orbit domain requires at least two, so **558** pairs pass that minimum-size gate. A residual orbit contributes exactly 166 subgroup points. An opposite 33 orientations have their small base entirely contained in the pool, but that fact alone does not measure relation yield. For example, two Gray-prefix K=600 streams can share 599 orbits, leaving just one residual orbit (166 points); hash and Gray-prefix streams in this archived panel are disjoint, leaving the pool's full orbit count as residual. Both extremes matter to the solver's domain size and partial-relation support, and neither by itself predicts runtime.

The archived panel is a finite set of **specific public inputs**, so the disjointness findings do not generalize to every seed, base size or hash policy. These local objects were downloaded and checked in the prior S3 replay cache; this stage made no S3 write or new live S3 claim. A representative K=600 Gray/hash cross-stream native import already passed separately, but its paired natural solver trial has not launched under the busy-host admission gate. This receipt measures no natural partial relation, useful rank, graph closure, linear algebra, individual logarithm or matched rho, and does not select a minimum-total-runtime factor base.

To reproduce, place the 27 committed S3 objects named `<compressed_blake3>.jsonl.gz` in one local directory, then run with Python `blake3` available:

```sh
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-pilot-overlap-20261015/overlap.py --objects-dir OBJECTS_DIR --check
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-pilot-overlap-20261015/audit.py --objects-dir OBJECTS_DIR --check
python3 research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-pilot-overlap-20261015/compare_prior.py --check
python3 -m unittest discover -s research/koblitz_n83_factor_base_sweep_20261008/verification/cross-stream-pilot-overlap-20261015 -p 'test_*.py' -v
```
