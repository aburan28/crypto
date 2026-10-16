# Native admission observations

Source `7e96187d3b819ecb6b5554d2024260bfe9eb0cba`; run started `2026-10-15T23:17:27Z`. Exact native counts from the pinned policy, plan, raw control entries and fresh certificates.

| Observation | Count |
|---|---:|
| Standard prime models | 54 |
| Construction adapters | 47 |
| Replay adapters | 47 |
| Candidate curve-degree pairs | 157 |
| Corrupted-artifact controls rejected | 8 |
| Resource stress accepted | 0 / 0 |
| Resource stress censored | 0 |
| Newly discovered maps | 0 |

P192/13 and P521/7 were admitted as separate known-answer controls. All eight altered artifacts and the forged saved-score submission were rejected. Every leaderboard build performed fresh replay.

[Native audit summary](audit.json), [complete audit](audit.full.json.gz), [policy](policy.json.gz), [candidate plan](plan.json.gz), [raw packet](raw.tar.gz), [diagram](admission.svg), [visual report](report.pdf). Compressed JSON files preserve their original bytes; the raw archive also contains every input, map, certificate and attempt receipt. Raw receipts are shared-host operational records; no speed comparison is made.
