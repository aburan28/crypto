# Independent n37 shared-rank replay result

The preregistered [Linux x86-64 CI run](https://github.com/aburan28/crypto/actions/runs/37178900016)
finished successfully. Its [full receipt](RECEIPT.json) has SHA-256
`bc2e2ebf42f248d55d81691ffb87d337eb0063c3c8e9bd516bd654bbde56729f`.
The auditor class `ECBENV2hda64b23436f4` differs from the source macOS arm64
class `ECBENV2hd82681268e96`; the auditor binary SHA-256 is
`e044499e84294b09d81b4e7459eda993732f59ceb142f8f76a6bfbb1767584c0`.
The receipt matches all five pinned source-session file hashes, passes every
integrity/regrading check, and reproduces **15/15 measured records** exactly;
all 18 executions in the session remain verified. The CI job itself passed.

The original measuring binary was not archived. The native `claim attach`
path therefore reconstructs the complete [original diagnostic claim](../CLAIM_FINAL_DIAGNOSTIC.json)
from the frozen session and requires exact JSON equality before adding any
independent metadata. It rejects a changed original online time and the
source host's same-class audit. The resulting [attached report](CLAIM_INDEPENDENT_DIAGNOSTIC.json)
has SHA-256 `adabcd19c30e81c6f3967377a13e8ef2acd23ade02eb26ca77b5fa2312d8ab5c`;
the native [schema check](CHECK.json) passes. Its two replay-certificate fields
name the same full receipt because that audit reproduced both IC and rho runs.
The pre-existing candidate, workload, run, target, phases, and timing fields
are unchanged. A structural JSON comparison finds exactly five changed keys:
`independent_validation`, `independent_replay`,
`independent_replay_pointer`, and the two arm replay certificates.

**Decision:** independent deterministic replay is satisfied for this one
session. The attached report remains an **L0 diagnostic**: its per-run numeric
`online_speedup` is a descriptive raw ratio, not an admitted performance
claim. The frozen [analysis](../RESULT_FINAL.json) keeps aggregate
`online_speedup: null`. The hosted VM audit did not repeat or isolate the
Mac wall timing, and both algorithms still leave native work unpriced. No
speedup claim or n131 transfer follows. The next measurement requires a
host-level isolation receipt (for example, isolab strict tier A plus ecbench
L3) for a same-point Linux session, common native-work pricing, and a
separately frozen new-target panel.

To regenerate the attached report, build the current source and run:

```sh
target/release/ecbench claim attach \
  --dir research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02 \
  --ic ic-shared --rho rho-strong --workload W157bda94ea05 \
  --base-report research/ecbench_n37_shared_rank_20261004/CLAIM_FINAL_DIAGNOSTIC.json \
  --independent-receipt research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/RECEIPT.json \
  --pointer research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/RECEIPT.json \
  --out /tmp/n37-independent-claim.json --exit-code
diff -u research/ecbench_n37_shared_rank_20261004/independent_validation_20261004/CLAIM_INDEPENDENT_DIAGNOSTIC.json /tmp/n37-independent-claim.json
```
