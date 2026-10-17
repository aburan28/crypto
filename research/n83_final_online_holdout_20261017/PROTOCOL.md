# Twelve sealed N83 points for the final online factor-base comparison

This fixture fixes twelve previously unused public subgroup points before
either pending Tier 1 factor-base panel has produced a selection. It reserves
an independent comparison of the smallest cold base and the strongest online
challengers. No timing, candidate choice, or resource result is inferred from
the point coordinates.

The exact curve is `EC1N83Ckb1a1haacb69ee9f37`: polynomial-basis
`GF(2^83)` with modulus integer `9671406556917033397649557`, equation
`y² + xy = x³ + x² + 1`, prime subgroup order `8569786107849059`, and
cofactor `1128547018`. The generator and subgroup checks are inherited
through the frozen fixture lineage. For index `i = 0..11`, the scalar law is
`1 + (little_endian(SHA256(domain || u32_le(i))) mod (r - 1))`, with domain
`n83-final-online-holdout-twelve-20261017-v1`. Only the public points are
producer inputs; [answers.jsonl](answers.jsonl) is an independent replay
sidecar. `generate.py --verify-sealed` regenerates the exact files, checks
subgroup membership, and proves non-overlap with the precision, finalist,
five-arm, K992 and earlier N83 fixtures represented by its parent chain.

| Frozen item | SHA-256 |
| --- | --- |
| [public targets](targets.jsonl) | `4d167f8b2b85bcf6c599494daf3b3a61210a9dc8c0d65ad3ebcdf07cd6d5c4c7` |
| [answer sidecar](answers.jsonl) | `f94a313ea6e7c414a4fc85c3ec55e384f6c6be14febbe3304d4b9dde1e54e484` |
| parent K992 generator | `514329656510e88d3316b688eef4d15104c64a5be7edb42f56ace269461d1dec` |

The final arm set is fixed by the following rule before any holdout timing:

1. Always include K992/root11
   `IC1N83Ckb1a1fb164672PDP4rootRCguidedLAgaussTDdirectISO0h829ed476f094`
   and K1107/root11
   `IC1N83Ckb1a1fb183762PDP4rootRCguidedLAgaussTDdirectISO0h82329aa062fa`.
2. Add the five-arm panel's exact online `selected_candidate_id` only if
   `online.status == SELECTED` and both host strata were admitted. Otherwise
   add its predeclared K1280/root11 control
   `IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dac`.
3. Add the K992/K1024/K1107 panel's exact online `selected_candidate_id`
   under the same admission rule. Otherwise add its K1024/root11 control
   `IC1N83Ckb1a1fb169984PDP4rootRCguidedLAgaussTDdirectISO0h3eb055926289`.
4. Deduplicate identical candidate IDs, preserving this order. All arms use
   the same frozen IC source SHA-256
   `8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`.

The later dispatch must bind those exact `IC1` manifests, the two panel
selection records, executed binaries, arm order, caps, Tier 1 resource
envelope, verifier and matched rho policy before any holdout measurement.
Run points 0–5 and 6–11 on independent Tier 1 instances, retain every
failed or noisy pair, and pool only independently admitted strata. The
primary comparison is each one-target IC online interval against an
independent rho solve of the same point. Select an exact online candidate
only if it is faster than every other arm on paired mean log time and every
corresponding exact two-sided `2^12` sign-flip test passes
`0.05 / binomial(number_of_distinct_arms, 2)`. Apply the same fixed rule to
paired cold-process millisecond differences as a separate supplementary
factor-base objective. An incomplete or noisy stratum leaves the selection
unknown. No point from this fixture is used to choose the candidate list or
tune an implementation before dispatch.

The implementation uses one fresh process per candidate and public point.
For point index `i`, rotate the deduplicated arm order by `i mod arm_count`.
Run an exposed-point rho calibration immediately before and after each IC/rho
pair; alternate IC then rho and rho then IC by parity of `i + fixed_arm_index`.
The paired rho seed is `2026102900 + arm_count*i + fixed_arm_index`.
The IC child cap is 480 seconds for `K < 2000` and 900 seconds otherwise;
each rho or calibration child has a 120-second cap. The producer cap is
190 minutes per host and the cloud job cap is 260 minutes. A host must pass
all 6 × `arm_count` same-point IC and rho arithmetic replays, their 12 ×
`arm_count` before/after calibration replays, exact source and base checks,
event order, resource telemetry and a ten-percent paired calibration spread
gate. The provider's distinct runner IDs, matching archived runner names and
overlapping job windows establish separate Tier 1 allocations. A shared
workflow job label is not used as a host identity.

`prepare_dispatch.py` writes `dispatch.json` only after both upstream panel
aggregates and their provider archive provenance have been versioned. A panel
that is not admitted contributes its prespecified fallback. The dispatch
file freezes the exact input hashes and arm choice. Only a push containing
`dispatch.json` to the dedicated
holdout branch starts the cloud workflow. Source and test commits before then
cannot start a measurement.
