# ECDLP harness run — 2026-09-29

Companion run of the Crypto Autoresearcher ECDLP harness against the
ECDLP-tagged eval problems and the core Semaev/rho experiment entry point.
This is an operational execution receipt, not a cryptanalytic advance.

## Scope

| Item | Value |
| --- | --- |
| Companion repo | `aburan28/crypto-autoresearcher` |
| Commit | `3e12d32e845a682d685db0a90c0a3847ea159058` |
| Host | Linux x86_64 cloud agent |
| Inference backends | unset (offline path only) |
| Unit | verified discrete logs / valid run records; no end-to-end IC `S` claim |

## Hypothesis (operational)

The Autoresearcher harness can:

1. materialize ECDLP eval fixtures and grade them offline with `harness.rho` /
   `harness.toycurve`;
2. complete `EXP-SEMAEV-001` (matched Pollard rho + S₃ Groebner measurement)
   with independent certificate verification;
3. recover a known 16-bit toy logarithm via DP rho and kangaroo walks.

Success: every selected ECDLP eval task passes its critical graders; every
Semaev/rho run record is `completed_valid`; walkviz recovers the planted `k`.

Stop: any critical grader failure, invalid certificate, or harness crash.

## Commands

```bash
git clone https://github.com/aburan28/crypto-autoresearcher.git
cd crypto-autoresearcher
pip install -e '.[dev]'
python3 -m orchestration doctor
python3 -m orchestration.eval validate --suite evals/suites/capability.yaml --backend local
python3 -m orchestration.eval validate --suite evals/suites/discipline.yaml --backend local
python3 -m harness.run --experiment EXP-SEMAEV-001 \
  --bits 8,10,12 --seeds 1,2 --factor-base 14 --out /tmp/EXP-SEMAEV-001
python3 -m harness.run_walkviz --seed 7 --field-bits 16 --out-dir /tmp/walkviz-b16
python3 -m pytest -q tests/test_harness.py
```

ECDLP eval problems were solved offline with `harness.rho` / exhaustive search
and graded through `orchestration.eval.graders.run_graders` (no model loop;
API keys were unavailable).

## Boundaries

- Floor: none claimed (toy-scale verification only).
- Reference: Pollard rho on the same toy instances (matched in Semaev runs).
- No ECC2K-130, m=83, or frozen IC tournament promotion claim.
- Agent-loop `orchestration.eval run` was not executed (no inference keys).
