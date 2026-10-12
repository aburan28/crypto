# Version 2: separate qualified cold and online references

This protocol continues the three-attempt campaign. Round one is sealed and is
replayed with its retained evaluator. Version 2 applies only to attempts two and
three, seeds 2026092552 and 2026092553. It does not reset the attempt count or the
familywise error budget. This contract update is not an improvement result.

## Accepted evidence and scope

The reference and observer evidence is the archive accepted in PR 876:
`ic-generic-reference-qualification-20260926.tar.zst`, SHA-256
`680369a9b822dcc3321380bff751bd2b94f0892c915f201a67916fc7bd49fc5a`.

The sorted-key compact JSON qualification object has SHA-256
`8479d9b9b777297d9c05dd860d5189ecc0791db9c70b6c641bc19ec986fae256`;
the observer object has SHA-256
`96b687f7887442f9059d733225274dabd07b25900263eabb17f721636089dee6`.
These are object identities, distinct from the byte hashes of the original files.

Use four frozen references. Names below are driver roles; exact per-cell IC1,
RHO1, curve, workload and run identities still come from admission.

| Role | Qualified selection | Source SHA-256 | Configuration |
| --- | --- | --- | --- |
| incumbent | incumbent | 240b8daa0478aacb1f6fc248de9fd5ed9b0773c59d3ae869ced0bab8be438378 | pair_table, tiny_gauss, m=3, batch_trials=1, max_trials=65536 |
| ic_online | prepared_both | 65df2c4417ec963a2b03e511719df28df51ea266a3c4a9e825307b8d1ac8c053 | same IC configuration |
| rho | rho_incumbent_8 | 240b8daa0478aacb1f6fc248de9fd5ed9b0773c59d3ae869ced0bab8be438378 | same configuration, rho_parallel_walks=8 |
| rho_online | rho_generic_dense_16 | 72da739cbc72e351e85c0184664ea60b6d2f23b54aab0e9ff45fdd1a60714d7d | pair_table, dense, m=3, batch_trials=1, max_trials=65536, rho_parallel_walks=16 |

The generic rho build identity is
`54316bcec91a68a5a08bcc8e92548a46ad75e2a27d8c056fe6520d79f3410962`.
Preserve its controlled build and executable, all dependency/source bindings,
and actual runtime dispatch. Prepared sources retain their source-bound
preparation and compiler checks.

All costs remain those of fully charged instrumented implementations. The
360-pair whole-mode observer study supports semantic agreement on its frozen
panel. It does not establish zero timer overhead or uninstrumented solver
performance. Legacy OBS1 observations remain unadmitted, with unknown scientific
phase costs. Do not subtract instrumentation or certificate costs.

The development selections supply fixed references, not evidence that one
prepared IC implementation significantly beat the other. Cold and online
leaders differ, so both remain references. Reference arms cannot enter the
challenger portfolio or win candidate selection.

## Workload and finite schedule

Use the original five development cells n17a1, n19a0, n23a0, n23a1, n31a0.
Confirmation adds n29a1, with twelve new public points per cell: 72 targets.
Each job has one supplied public point and three process repetitions. Replay
uses the same 72 confirmation points in fresh processes. Never tune using
confirmation or replay data.

Before preparing a round, verify every prior round with its own frozen evaluator
and carry its exact contract/decision chain. Extend the original 5,133-point
history with every exposed point from preceding rounds and supplemental
qualification, readiness and control panels, including all 25 points in the
new reference study. An altered seed cannot make an old public point fresh.
Keep failed and interrupted preparation evidence in the exclusions.

A candidate panel must be registered before execution, with exact configurations,
source/build identities, hypotheses and stopping conditions. Version 2 permits
at most eleven IC candidate arms including the incumbent: ten challengers.
This is below the original maximum of sixteen so that the additional reference
fits the unchanged resource budget. Retain six diverse challengers when eligible,
including an exploration slot; preserve online/cold leaders, families and
measured combinations of mechanisms. Failed candidates remain recorded.

At the maximum panel size, the schedule is:

| Stage | Public points | Arms per point | Repetitions | Native/profile pairs |
| --- | ---: | ---: | ---: | ---: |
| A/A | 5 | 2 | 3 | 30 |
| Smoke | 5 | 14 | 3 | 210 |
| Development | 15 | 14 | 3 | 630 |
| Selection | 15 | 10 | 3 | 450 |
| Confirmation | 72 | 5 | 3 | 1080 |
| Replay | 72 | 5 | 3 | 1080 |
| Total | — | — | — | 3480 |

A/A compares the cold incumbent with a copy. Later stages always retain both IC
references and both rho references. Smaller or failed portfolios reduce the
executed schedule; they never create successful-subset estimates or permission
to add replacements after viewing results. The cap stays 3,500 native/profile
pairs, Linux amd64, Rust 1.94.1, Valgrind 3.22.0, one pinned CPU and one Rayon
thread, 8 GiB and 180 seconds per child. Failures, timeouts and OOMs remain rows.

## Accounting, inference and decision

The primary metric remains native online wall time for one previously unseen
supplied public point, after reusable preparation through scalar replay.
Preserve all unsuccessful target-dependent work. Complete cold native process
time and complete user-space guest instructions are supplementary acceptance
metrics. Every phase stays exclusive and fully charged; missing totals remain
unknown. Report S=Ir/sqrt(r), actual usable base count B, folded columns K,
paired rho ratios and the existing weak implementation-specific K-instruction
collector floor where applicable. The floor is not a universal IC lower bound.

Pair cold instructions and cold process time with `incumbent`; pair online
time with `ic_online`. Development retention and selection use these explicit
metric baselines. Online IC/rho tables still compare each candidate with the
same-point rho observations directly.

For confirmation and replay, use the original fixed-cell paired-target basic
bootstrap, 50,000 draws, seed 914223 + attempt number. Reuse the same sampled
target indices across the three metrics, including when their baselines differ.
The family remains three attempts × two final stages × three metrics = eighteen
tests, with one-sided alpha 1/360 and conservative Monte Carlo tail index 137.
Coverage remains nominal and approximate, conditional on the declared target law
and fixed cells. A process repetition is not a new target.

Promote only if both final stages independently verify all required targets and
references, cold instruction and native-time ratios are each at most 0.8, online
ratio is at most 1, every one-sided familywise upper bound is below 1, and every
per-cell ratio is at most 1.1. Preserve the six-cell allocation and all eighteen
tests; changing references does not grant a new inference budget. If a reference
is incomplete, the comparison is inconclusive.

Deliver the immutable run record, sources, profiles, native outputs,
certificates, exact commands, failures and independently replayed decision in
an evidence PR, with the same results in the canonical scoreboard. No broader
curve-family, asymptotic or global-optimum claim follows from this finite panel.
