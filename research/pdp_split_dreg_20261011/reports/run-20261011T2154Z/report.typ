#set page(paper: "a4", margin: 18mm)
#set text(size: 9.5pt)
#set heading(numbering: none)

= Run report: matched-rho ecbench session `rho-01`

*Study:* `research/pdp_split_dreg_20261011` (PROTOCOL.md, "Matched rho"). *Descriptive only:* this page records what ran and what was measured; interpretation is in the study README.

== What ran

- *Spec* `SPEC-rho.json`, id `ECS1h10c3ed7d3628`; *session* `ECBS1he23809badf5b`, directory `sessions/rho-01/`.
- *Command:* `ecbench run --wait --spec SPEC-rho.json --out sessions/rho-01 --cpus auto`.
- *Commit* `51bb1dc5e1fea5a7b53a2c6c8b02f51f75e77233` (clean); *binary* SHA-256 `909efb98…3d5740`.
- *Host:* Intel Xeon @ 2.10 GHz, 4 logical CPUs, 1 NUMA node, virtualised Linux 6.18 (Claude Code cloud container); runner CPUs 0–2. Started 2026-10-11 21:54:49 UTC, finished 12.4 s later.
- *Arms:* `rho-strong` = `rho.signed_frobenius_strong` (lanes 32, dp_bits 4, step_cap_factor 2000; reference); `rho-frob` = `rho.signed_frobenius` (cap_multiple 64; baseline).
- *Workloads:* Koblitz `icv1-f2m11-t67-05f5aa36`, `icv1-f2m13-t181-515ee569`, `icv1-f2m17-tm101-00378d4e`, `icv1-f2m19-tm797-9c54981b`; 8 public one-target workloads each (target seed 202610111101); 5 measured rounds + 1 warm-up; alternate order, seed 202610111102.

== Completion and isolation

384 executions planned, 384 recorded, all `verified` (status counts `{"verified": 384}`). Isolation earned: L2 on most runs, L1 on 3–11 per (arm, curve) cell (blockers to L3: CPU not isolated, no performance governor). Wall time is not used by this study.

*Audit:* `ecbench verify --replay-all` → `audit OK: 384 records, 384 verified, 320 replays reproduced, 0 problems`; receipt `audits/rho-01.audit.json`, SHA-256 `0b56a9b730d9d14b726392d15500ffd22093413f048a72c4d1171099c33dcc0c`.

== Measured values (ecbench table, unit S = charged group-addition equivalents / √r; all lower bounds by ecbench's rule, canonicalisations uncharged)

#table(
  columns: 6,
  align: (left, left, right, right, left, right),
  [arm], [curve], [log₂ r], [mean S], [95 % interval], [S / floor],
  [rho-strong], [f2m11], [9.95], [42.024], [21.121–77.100], [157.3],
  [rho-strong], [f2m13], [10.97], [16.650], [14.037–23.827], [67.7],
  [rho-strong], [f2m17], [16.00], [4.054], [3.913–4.201], [18.9],
  [rho-strong], [f2m19], [18.00], [2.378], [2.289–2.460], [11.7],
  [rho-frob], [f2m11], [9.95], [2.448], [2.366–2.575], [9.16],
  [rho-frob], [f2m13], [10.97], [2.095], [1.752–2.688], [8.52],
  [rho-frob], [f2m17], [16.00], [1.093], [0.914–1.330], [5.09],
  [rho-frob], [f2m19], [18.00], [0.577], [0.515–0.661], [2.84],
)
n = 40 verified runs per row (8 workloads × 5 rounds). Source: `audits/rho-01.table.md`; per-workload pairs in the session's saved comparison.

#figure(image("rho-session-s.svg", width: 100%), caption: [Mean S with ecbench's 95 % interval per arm and curve, session `ECBS1he23809badf5b` (recorded values; source `audits/rho-01.table.md`, layout `reports/plot_rho_session.py`).])

== Raw records

`sessions/rho-01/records.jsonl` (384 lines), `sessions/rho-01/session.json`, `spec.json`, `plan.json`, `host.json`; audit and table under `audits/`.
