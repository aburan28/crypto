# Volcano experiment downloads, 17 October 2026

[Open the Drive folder](https://drive.google.com/drive/folders/109DpcQpX7ATz6GeI8iSWAExGMWBoR-H_).
These are frozen snapshots. Complete graph enumeration and completed ECDLP
measurements have different coverage.

| Artifact | What it contains |
| --- | --- |
| [Interactive graph explorer](https://drive.google.com/file/d/13AQHFuLnePIfbbS9XbaXNueJY-JMyKCE/view) | All 24,823 vertices, 51 conductor strata and 69,300 recorded directed kernel maps. Download the HTML and open it; Drive preview does not run it. |
| [99-vertex measured archive](https://drive.google.com/file/d/1hLVWVad2SrgHaBn6U5iSyqnKbx5cAkeY/view) | Complete three graph databases, exact frozen dev binaries/source, full 24,823-row queue, execution ledger, 99 sealed sessions, replay receipts, native statistics and report. |
| [Three-page report](https://drive.google.com/file/d/1h_G5uJZYYqVzBL3_kynRdYZvE8vjJu9z/view) | Measured coverage and yield/cost summaries over all 51 strata; every page inspected. |
| [Detailed statistics CSV](https://drive.google.com/file/d/1OWEkYFgWuK_wgMu0auEHsYxJ2VDugFIA/view) | Per-graph/stratum/arm/outcome min, max, arithmetic mean and sample count, including phase/native-operation counters. |
| [Native pricing band](https://drive.google.com/file/d/1NV-4hpuoG23Xf0TitT8yUBPGAnRJCoj2/view) | All 792 curve/arm rows at common unpriced-operation weights 0/0.1/1 GAE; incomplete/unaccounted cases remain explicit. |
| [Separate optimized control/pilot](https://drive.google.com/file/d/1lJt6NhNwXVpvdSBoXqCXhoi-flqEzwmF/view) | Actual release runs, binaries, exact dev baseline capsules and native comparison receipts. This is not the optimized whole-queue survey. |

The three source models are `icv1-fp57-t536870914-6cabb0c9`,
`icv1-fp57-tm536870941-e6d0c766` and `icv1-fp59-t1073741954-4dd0b1ef`.
Their maximum endomorphism conductors are 8192, 6561 and 1728, respectively.
Their largest prime subgroups have orders 61751, 19768889568671 and 4253;
field size is not subgroup security size. Every vertex uses its certified
actual largest-prime subgroup, not a substituted smaller subgroup.

The stopped dev checkpoint has **99 measured vertices across all 51 strata**:
1,584 attempts, 1,065 verified, 422 exhausted and 97 timeouts. It compares
projected factor bases 16/64, MITM/subtraction, a raw-base control and matched
rho/BSGS, using two public targets, one round and bounded trials/watchdogs.
Every completed successful/exhausted deterministic report has a replay receipt.
These are same-implementation replays, not independent validation.

It found 57,176 relations from 2,223,204 source trials and 2,075,145 actual
decomposition calls. Repeated queries skipped before the oracle explain the
different denominators; summary v2 exposes both. Timeout counters are unknown,
not zero. The old first-12 report's "decompositions" label meant source trials;
historical records/reports/uploads remain unchanged. Do not pool different
subgroups, build profiles or repeated control/pilot targets as one population.

The optimized control/pilot have 60 attempts: 40 verified, 18 exhausted and
two timeouts. All 56 comparable completed deterministic profile pairs match
exactly; four baseline-censored pairs have no final counts to compare.
The optimized full-queue survey has its own binary/ledger and 24-hour/2-GiB
guard; its execution status is in `execution/optimized-continuation/` and the
local ledger. Its results are separate from this dev snapshot. No runtime
speedup, fastest-method or asymptotic claim is established by these runs.

The graph CLI is `crypto isogeny class`. For a supported ordinary prime-field
curve JSON, use `plan --curve curve.json`, then
`build --curve curve.json --out graph.sqlite --steps all`. Inspect with
`vertex --input graph.sqlite --id N` and `status --input graph.sqlite`.
The native [`CLI documentation`](https://github.com/aburan28/crypto/blob/codex/full-isogeny-class-cli-20261016/tools/volcano-cli/README.md)
defines backend limits, resume, map verification and the separate one-target
BSGS command. Graph generation does not itself run graph-wide index calculus.

Graph completeness is for Fp-isomorphism vertices and the declared finite
degree sets, not other primes or infinitely many compositions. The original
A–E graphs/campaign, full eight-target/five-round/million-trial matrix, exact
per-vertex representation registration, transported targets, composite PH and
independent replay remain open. Finite-field IC for the actual embedding
degrees is not implemented; the k=2 MOV teaching example is not that backend.

Local extraction/hash checks and deterministic archive repacking pass. Drive
metadata confirms names, sizes and parent folder; remote whole-file hashes
were not checked. Delivery JSON files retain exact byte counts and SHA-256s.
