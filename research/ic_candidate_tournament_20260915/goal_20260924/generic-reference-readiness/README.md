# Generic readiness on the five development cells

All fifteen registered jobs completed and passed independent admission: ten IC
solves (dense and sparse relation-LA policies on each cell) and five rho solves.
This establishes admission readiness for these pair-table configurations. It is
not reference qualification, a speed comparison, or an improvement round. The
incumbent is unchanged; one round is closed without a qualifying winner and two
remain. The implementation depends on PR 822, which was unmerged at registration.

The [protocol](PROTOCOL.md) was committed as
`9faa3a2ce2651cee6a9ac3cf94d179cd6cc28d81` before execution. It freezes one public
target per cell, seed 2026092561, two IC LA policies, requested rho width four,
60-second child deadlines, one Rayon worker and a 20-minute driver deadline.
The measured producer is the unchanged source-bound final build from generic
scientific admission. The local host is macOS arm64; unavailable affinity and
address-space controls remain null. The outer driver exited successfully after
11.899312 seconds; this includes external auditing and is not a solve-time result.

## Stored result and interpretation

The [result records](RESULTS.json) retain every canonical IC candidate/workload/run
key. Each IC row independently reconstructs the actual usable base, stored
matrix and full rank, query law, descent witness, scalar replay and phase closure.
The same supplied point is used by all three pipelines in each cell.
The [curve-identity audit](curve-identity-audit.json) reads the hash-verified
qualified-reference archive and confirms all five exact canonical curve IDs,
including field representation, subgroup and generator. Matching degree labels
alone would not establish that compatibility.

| Cell | Actual usable base points | Folded columns / final rank | Collection attempts, dense / sparse | IC completions | Rho completion |
|---|---:|---:|---:|---:|---:|
| n17a1 | 272 | 8 / 8 | 8 / 8 | 2 / 2 | 1 / 1 |
| n19a0 | 304 | 8 / 8 | 8 / 8 | 2 / 2 | 1 / 1 |
| n23a0 | 368 | 8 / 8 | 8 / 8 | 2 / 2 | 1 / 1 |
| n23a1 | 368 | 8 / 8 | 12 / 12 | 2 / 2 | 1 / 1 |
| n31a0 | 496 | 8 / 8 | 10 / 10 | 2 / 2 | 1 / 1 |

Both n23a1 IC runs preserve four **unresolved** collection attempts. These are
algorithm outcomes, not mathematical proofs that decomposition is impossible.
There are no `proved_unsat` claims in this panel. The checker's existing bound
on independently proved negative PDP results therefore remains unchanged; these
controls do not qualify larger negative proofs for F4/F5/SAT or other engines.
Each target descent has one witnessed attempt. Sparse policy admission does
not imply that its filtered small matrix executed a nonempty Wiedemann core.

All ten IC run keys are unique. All eleven existing environment-override controls
were rejected. No failure was retried and no worker was rerun during evidence
replay. Performance qualification, online speedup, normalized S and promotion
remain unset. There is no rate estimate or timing ranking from five correctness
fixtures. Calibrated native/instruction comparisons with strong references and
instrumentation controls remain the next gate.

## Reproduction and future exclusions

[EVIDENCE.json](EVIDENCE.json) binds a 1,734,149-byte archive containing 171 files
(4,798,195 expanded bytes), including every raw report, canonical record, worker,
build receipt, source manifest, frozen Python checkers and driver log/exit.
Archive SHA-256:
`6eb4327985705336aadbeb7c981010c4c0469655fccc384f1e1972c570236cb0`.
The unchanged worker's Rust source archive is identified in that manifest.
[Fresh extraction](fresh-restore.json) verifies every file hash and reproduces
the result JSON byte-for-byte, without executing the archived worker.
The full [Python suite](tests.log) runs 166 tests (164 pass, two existing local
platform skips); the [site suite](site-tests.log) passes 39 tests. The added
regression extracts this archive, verifies its full manifest and runs its frozen
checker in a separate process to reproduce every result record.

```sh
mkdir /tmp/ic-readiness-restore
tar -xzf research/ic_candidate_tournament_20260915/evidence/ic-generic-reference-readiness-20260925.tar.gz \
  -C /tmp/ic-readiness-restore
python3.12 /tmp/ic-readiness-restore/ic-generic-reference-readiness/checkers/goal_20260924/generic-reference-readiness/replay.py \
  --evidence /tmp/ic-readiness-restore/ic-generic-reference-readiness \
  --out /tmp/ic-readiness-replayed.json
```

The five public points are retained in [fixtures.json](fixtures.json). Include
all five in the next registered target-exclusion history, along with every point
from the prior rounds and other development studies. The existing `target_history.extend`
can consume this directory. Do not overwrite the sealed round-one history or
reuse these exposed points as fresh confirmation targets. The held-out n29a1
cell was not used.
