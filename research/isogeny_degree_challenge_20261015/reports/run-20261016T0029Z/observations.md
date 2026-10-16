# Native admission observations

Source `9bb7cbce2b84151eff6846b02da842a8ef019343`; run started `2026-10-16T00:29:57Z`. Exact native counts from the pinned policy, plan, raw control entries and fresh certificates.

| Observation | Count |
|---|---:|
| Standard prime models | 54 |
| Construction adapters | 47 |
| Replay adapters | 47 |
| Candidate curve-degree pairs | 157 |
| Corrupted-artifact controls rejected | 9 |
| Resource stress accepted | 20 / 20 |
| Resource stress censored | 0 |
| Newly discovered maps | 0 |

P192/13 and P521/7 were admitted as separate known-answer controls. All listed altered artifacts and the forged saved-score submission were rejected. Every leaderboard build performed fresh replay.

The starting leaderboard re-admitted these historical maps. All passed twenty public scalar transports per map.

| Exact source | Degree | Exact target | Coverage |
|---|---:|---|---|
| `icv1-fp192-t31607402316713927207482677199-52e4af59` | 10453 | `icv1-fp192-t31607402316713927207482677199-e2f29854` | PARTIAL |
| `icv1-fp224-t4733100108545601916421827343930821-fec01e99` | 1471 | `icv1-fp224-t4733100108545601916421827343930821-5a35f6f7` | PARTIAL |
| `icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-f767e0d1` | 7 | `icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-5014bab6` | COMPLETE |
| `icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-f767e0d1` | 7 | `icv1-fp521-t657877501894328237357444332315020117536923257219387276263472201219398408051703-e804740d` | COMPLETE |

The one-second P192/10453 control returned `CENSORED` / `TIMEOUT` with no score. The full-budget leaderboard used a separate fresh replay.

[Native audit summary](audit.json), [complete audit](audit.full.json.gz), [policy](policy.json.gz), [candidate plan](plan.json.gz), [raw packet](raw.tar.gz), [diagram](admission.svg), [visual report](report.pdf). Compressed JSON files preserve their original bytes; the raw archive also contains every input, map, certificate and attempt receipt. Raw receipts are shared-host operational records; no speed comparison is made.
