# Third registered generic-backend campaign: negative F4/F5 smoke qualification

The single authorized `workflow_dispatch` is [Actions run 37399021219](https://github.com/aburan28/crypto/actions/runs/37399021219), attempt 1, on main checkout `63a09ff7b297531b8bacef853d79745b3575983b`. The frozen worker is `765c3c5f19032bd852163805f257c56babef2040`. Panel seed `2026093001`, panel SHA-256 `df92d5507785446a2a5b333bd7a04776781de4f54c63c5ea3ffa921bc99ad2e4`. Preflight and artifact smoke passed. The campaign job finished, and `tournament.py verify` inside that job reported 45 receipts verified. Seed `2026093001` is consumed.

The packed archive is GitHub artifact `ic-generic-backend-qualification-v3-37399021219-1` (artifact id 11387084743, 30-day retention). Its zstd payload SHA-256 is `70e9de519a36742f1a2d47945b5dfbda13f4e43c0f4886f36e9f0e748484f851`, 607,892,011 bytes, matching the pack manifest and a local rehash. Uncompressed tree size reported by the packer is 2.19 GiB. This note retains the gate, summary, natural-yield audit, and the ten generated fixtures extracted from that archive. It does not replace an independent replay of the full tree.

| Retained file | SHA-256 |
| --- | --- |
| `retained/family-gate.json` | `ed055283d42575a2e672e5eef20fea942002958c72405e36427ee65712b76ee2` |
| `retained/summary.json` | `75571e24d860ff59af3e37bec0051fb8731447def851c97f5f4df7c734d8fc94` |
| `retained/natural-yield.json` | `6d476045189f0fc4717e68781e0f84dc035d5a4d7af452a6e61475719b550264` |
| `retained/fixtures.json` | `66206b7fb050c4cda8709534490ca1b6b98abf0b53624c1e1353a9b6e2299eac` |
| `retained/archive-manifest.json` | `ab5f085304b5d3322423d1f346a056679e6260398427a59ea172cd11fe0b47ed` |

## Gate

`generic_backend_gate_v3.py` returned `NEGATIVE_FAMILY_QUALIFICATION`. `qualified_f4_f5_arms` is empty. Summary status is `AUDITED_WITH_FAILURES`: 45 trial slots, 34 verified native/profile pairs, 11 retained failures, 10 exposed points, `promotion_eligible` false. SAT is out of scope. This is not an IC-versus-rho speedup, a held-out confirmation, or a family promotion.

| Arm | Smoke scheduled | Verified complete | Censored | Ordinary queries retained |
| --- | ---: | ---: | ---: | ---: |
| `generic_f4_subspace_dense` | 5 | 0 | 5 | 0 |
| `generic_f5_subspace_dense` | 5 | 0 | 5 | 0 |
| `generic_pair_subspace_dense` | 5 | 4 | 1 (`n31a0`) | 14,847 on the four complete processes |

Every F4 and F5 cell has yield qualification `incomplete`. The natural-yield auditor records the censor reason `invalid or interrupted profiler JSON`, zero observed ordinary queries, and unknown rate (`no ordinary query exposure; rate is unknown`). The job progress lines marked those same processes `TIMEOUT` with null operation counts. That is a censored process, not a measured zero-yield solve and not an UNSAT proof. Neither algebraic arm can be qualified by dropping those jobs.

The pair-table arm is the subspace continuity control. It completed four cells and was censored on `n31a0`. Its query counts and instruction totals are not an F4/F5 result and are not a speedup claim. A/A `incumbent` and `aa_control` verified on all five cells. Smoke `incumbent`, `ic_online`, `rho`, and `rho_online` also verified on all five cells.

The ten exposed public points, five A/A and five smoke, are in `retained/fixtures.json`. Any later panel must exclude them, together with the two earlier censored corpora and the three sealed improvement rounds. Do not redispatch seeds `2026092901`, `2026092902`, or `2026093001`.
