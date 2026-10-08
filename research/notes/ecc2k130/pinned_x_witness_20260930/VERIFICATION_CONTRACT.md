# Separate the selected point lift from the pinned x-only witness

The compact S3 query reports two kinds of evidence for a target `Q`. The
selected four factor-base points prove a group sum and a recovered logarithm.
The two `pinned_intermediates` are x-coordinates of pair sums found earlier by
the x-only root search. These witnesses need not use the same signs when two
factor-base points share an x-coordinate. The frozen 2026-09-30 verifier
required them to be the same lift and consequently rejected one valid public
Q at n37/L1024 in every compact arm. Its failure and all six timing
ineligibility decisions remain unchanged.

For a future protocol, [check_witness.py](check_witness.py) checks both claims
independently:

1. Each selected index identifies a curve point with the reported x-code and
   a checked representative label. The four selected points sum to `Q`; their
   logged coefficients sum to the recovered scalar, which equals the
   separately verified fixture scalar.
2. Among the 16 sign lifts of those same four x-codes, at least one pair of
   pair sums has the two reported pinned x-coordinates as an **unordered
   multiset** and its four-point sum is `Q`. On this binary curve, the only two
   points with a nonzero x-coordinate are `P` and `-P`, so the search is exact.
   Infinity pair sums are rejected. The returned sign bits make the witness
   reproducible.

The [archived regression replay](replay_archived_failure.py) rehashes the
original n37/L1024 compressed artifact and 61 selected raw members against
the committed manifest. It independently replays all 15 full-rank traces
and checks fixture scalars, then
replays target indices 0 and 763 in each of 15 compact arms. All 15 instances
of Q763 have a valid pinned x-only witness but a different selected pair-root
multiset; ordinary Q0 matches the selected lift. The prior verifier still
rejects Q763. Mutating an in-field pinned root, out-of-field pin, recovered
scalar, or selected point index is rejected by the new check. The deterministic
[receipt](REPLAY.json) is a correctness regression only, not a new experiment
or timing admission.

Run from the repository root:

```sh
python3 research/notes/ecc2k130/pinned_x_witness_20260930/replay_archived_failure.py \
  --out /tmp/pinned-x-witness-replay.json
cmp /tmp/pinned-x-witness-replay.json \
  research/notes/ecc2k130/pinned_x_witness_20260930/REPLAY.json
```

The next cold protocol must freeze a new disjoint public-Q corpus, source and
input hashes, this witness contract, and an auditable Linux isolation gate
**before** timing. It must report the one-target online intervals required by
`AGENTS.md` separately from a batch-throughput question. The old hosted run
cannot be promoted by replaying it with this code. Its 0/6 timing-admission
result and the scoreboard's unset method-level speedup stand.
