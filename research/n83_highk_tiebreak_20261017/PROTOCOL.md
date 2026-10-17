# N83 high-K factor-base tie-break, frozen before either comparison result

The current IC source has two independently replayed high-K challengers:
K3450/root11 wide3 Bloom,
`IC1N83Ckb1a1fb572700PDP4rootRCguidedLAgaussTDdirectISO0hbde4a785a36e`,
and K3500/root4 standard Bloom,
`IC1N83Ckb1a1fb581000PDP4rootRCguidedLAgaussTDdirectISO0hb3a616c4e22f`.
Their exact manifests, source/binary hashes, full-base arithmetic replays and
exposed-point cost controls are in
`research/n83_k3072_current_source_20261017/RESULT.md`. This comparison
tests which base gives the shortest verified **one-target online** interval
under the same source and host envelope. The complete cold process is a
supplementary objective.

The third arm is the **online** selection of the separately frozen
K3072/K3450-dynamic/finalist ten-point comparison, if both of its Tier 1
strata are admitted and its online selection status is `SELECTED`. If that
gate is null, the fixed anchor is current-source K3072/root11 dynamic,
`IC1N83Ckb1a1fb509952PDP4rootRCguidedLAgaussTDdirectISO0hec7ae07d8b9a`.
The fallback is a comparator, not an inferred winner. Freeze and
independently replay the predecessor's validated rows, selection, and bound
provider archives before dispatch. The eventual dispatch must bind their
SHA-256 digests and the three exact candidate IDs. No outcomes from this
fixture may change the arm rule.

Ten public N83 subgroup points were fixed by
`SHA256("n83-highk-tiebreak-ten-20261017-v1" || u32_le(index))` for indices
0–9. Their [public-point file](targets.jsonl) has SHA-256
`c31ed8c9b8ceec71cf4570919bf01aa9bb8168077e37851a9ec144a55058a4b2`;
the [answer sidecar](answers.jsonl) has SHA-256
`5c104cfd26be3077b3b0081e46547d2b4c7a040809bd61b2e936198774422fd6`.
The sidecar is for independent replay and is not a producer input. The
[generator](generate.py) verifies the earlier fixture chain, point
uniqueness, curve membership, subgroup membership, and known-scalar
construction. The [read-only inventory](scan_prior.py) of 518 prior N83
JSONL files found no overlap among 32,376 earlier two-string point rows.
The ten points are
disjoint from the finalist and K3072 followup fixtures and stay unopened by
the producer until this protocol and the runner are committed.

Use two dedicated Tier 1 cloud allocations, each with five points. On each
point, run all three IC arms in fresh processes with empty base, index, and
log caches. Pair **each** IC arm with an independently verified strong rho
solve of the identical public point. Rotate the three IC arms by global
point index. Alternate `IC,rho` and `rho,IC` by global point and fixed arm
index. Run the frozen exposed-point rho calibration before and after each
pair, using seed `2026101601`. Set each paired rho seed to
`2026102100 + 3*point_index + fixed_arm_index`; indices 0, 1, and 2 mean
wide3, root4, and anchor, respectively. Start the IC online clock at the
first target-dependent computation after precomputation and stop after
independent scalar verification. Start and stop rho at the analogous
target-dependent boundaries. Keep the five exclusive IC online phases,
full cold process, all failed attempts, memory, and source/target receipts.

Use a 1,200-second cap for each high-K IC process, a 1,200-second cap for
the anchor if it is K3450 dynamic, 900 seconds if it is K3072, and 480
seconds for a finalist anchor. Cap rho and calibration processes at 120
seconds each. Cap each five-point producer at 240 minutes and each cloud
job at 290 minutes. Preserve every timeout, failed replay, noisy pair,
resource event, and global-cap stop. Do not classify a timed-out arm as a
slow but successful result.

Each host must independently pass the Tier 1 gate in
`docs/IC_BENCHMARK_TIERS.md`: immutable allocation/source/input and resource
envelope; exposed topology, affinity, cgroup memory and throttling; zero
observed steal time for admitted pairs; calibrated host noise; exact event
order; complete same-point arithmetic and scalar replay; and no missing
pair. Each before/after rho calibration must differ by at most 10% of its
mean. Pool strata only after both pass. Retain all raw archives and bind
their provider digests to independently replayed local copies.

For each candidate, report ten paired `rho_online_ms / IC_online_ms` values,
the median and each host's mean, and the paired 95% interval for rho minus IC
online time. A Tier 1 speedup on a host requires IC online below rho on all
five of its points. Compare the three IC arms on all ten matched points by
the mean log ratio of their online intervals, descriptive paired 95%
intervals, and exact two-sided sign-flip p-values over all `2^10` sign
assignments. Select an online winner only if it is faster in the observed
direction than both other arms and each of its two p-values is at most
`0.05/3`, accounting for all three unordered pairs. Otherwise record
`BOUNDED_UNKNOWN`. Apply the same three-pair gate to complete cold-process
costs and report that selection separately. This staged comparison tests
the high-K challengers against the previous selected leader under an
identical point and host envelope; it does not infer cross-host timing from
separate panels.
