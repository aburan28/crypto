# N83 root11 factor-base neighborhood screen

This bounded screen tests whether a nearby orbit count improves the selected
K3300/root11 factor base under the same exact N83 four-summand root index.
The ordered arms are **K3200, K3250, K3300, K3350**, each using five-byte
`2^30` root slots, reconstructed 11-byte states, rank seed 29, four rank
threads, one target thread, a fresh process, a 14 GiB RSS guard and a
1,200-second child cap. The source and fixture hashes are checked before
the first arm. An arm that reaches a cap, fails replay or exceeds the guard
remains a recorded status and cannot win selection.

All arms use the same 4,000 public subgroup points, in the same order, from
`research/n83_factor_base_confirm_20261015/targets.jsonl` (SHA-256
`451df8524b2c356dd4a523260620c4dff406975a122762b98ddddb7987b75c2c`).
The scalar answer sidecar SHA-256 is
`8e18c9c97ba9f57e5eb43f1b39b9e9575b5cfab87840c25c428dac150c868101`.
These points were previously used to compare K3300/root11 against
K3100/root22. This is a **secondary exploratory neighborhood screen**;
any newly selected base requires fresh one-target Tier 1 IC/rho validation.

The curve is `EC1N83Ckb1a1haacb69ee9f37`, polynomial-basis `GF(2^83)`
with defining polynomial bitmask `9671406556917033397649557`, curve model
`y²+xy=x³+x²+1`, subgroup order `8569786107849059`, and cofactor
`1128547018`. Each exact base is the first K distinct signed-Frobenius
subgroup orbits in the frozen construction; its actual usable point count
must be `166K` and its final log rank must be K. The K3300 source base hash
must equal
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`.
The producer source SHA-256 is frozen as
`a1eb3b937106bb7cfdee6c54d6422ceb126b628462361c9ff7a127127f230c89`.
The shared parent runner still carries the prior source pin; this study
explicitly checks that old pin and replaces it with the new exact SHA-256
before invoking the parent's input verifier. The intervening source edit
boxes wide-field state, and the K3300 base replay remains an exact check
against the prior N83 base digest.

For each K other than 3300, compute the paired differences
`online_ms(K, i) − online_ms(3300, i)` over all 4,000 points. Partition the
ordered target list into 80 contiguous blocks of 50, resample blocks
10,000 times with seed `2026101599 + K`, and take the sorted 2.5% and
97.5% endpoints. An arm is confirmed faster only when its upper endpoint
is below zero. Select the confirmed arm with the smallest measured mean;
otherwise retain K3300. Reproduce the exact result on a new target panel
before promoting a changed factor base to a primary Tier 1 comparison.

The online interval includes target query, point decomposition, relation
check, descent and recovery check after target-independent base/index/log
preparation. Record setup, final rank, root-key count, table/state bytes,
full process wall, peak RSS, status for every arm and every target's phase
readings. The independent verifier replays the full base and every
four-point relation and scalar before applying the frozen selection rule.
Hardware measurements are labeled `cloud_instance_exploratory`; they do
not inherit a Tier 1 wall-time claim from the earlier single-target run.
