# ECDLP harness runs — 2026-10-01

Seven consecutive executions of the Crypto Autoresearcher ECDLP harness
offline protocol (same commands as
`research/ecdlp_harness_run_20260929/PROTOCOL.md`). Operational receipt only;
no cryptanalytic advance, no algorithm change in this repository.

## Frozen commands

```bash
cd /Volumes/SSD990/crypto-autoresearcher
python3 -m orchestration doctor
python3 -m orchestration.eval validate --suite evals/suites/capability.yaml --backend local
python3 -m orchestration.eval validate --suite evals/suites/discipline.yaml --backend local
python3 -m harness.run --experiment EXP-SEMAEV-001 \
  --bits 8,10,12 --seeds 1,2 --factor-base 14 --out <OUT>
python3 -m harness.run_walkviz --seed 7 --field-bits 16 --out-dir <WALKOUT>
python3 -m pytest -q tests/test_harness.py
```

Environment: macOS (darwin) arm64, Python 3.12.8, offline path (no model
loop, no API keys). Out dirs were fresh per execution under
`/Volumes/SSD990/llm/tmp/opencode/` (the harness refuses to overwrite
immutable run records, so reruns used `-rerun`…`-rerun6` suffixes).

## Success / stop conditions

Success: `orchestration doctor` Ready; both suites validate; every
EXP-SEMAEV-001 run record is `completed_valid` with a verified certificate;
walkviz trace closes; `tests/test_harness.py` fully passes.
Stop: any critical failure, invalid certificate, or harness crash.
