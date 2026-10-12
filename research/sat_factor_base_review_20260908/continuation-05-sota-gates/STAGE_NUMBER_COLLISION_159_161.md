# Stage-number collision preservation

The `n=59` collector chain on `main` and the native-F4 branch independently assigned Stage numbers 159–161 before their histories were merged in PR #660.

The canonical native-F4 Stage 162–170 audit chain pins its own historical Stage-159 audit at `stage-159-current-gate-audit-20260922` and the native-F4 composer scripts `compose_koblitz_stage159_gate_audit.py` through `compose_koblitz_stage161_gate_audit.py`. Those pinned paths therefore retain the native-F4 artifacts.

The displaced `main` collector evidence is preserved additively:

- `stage-159-bucket25-current-gate-audit-20260922/` contains the original collector Stage-159 audit and seal.
- `scripts/compose_koblitz_stage159_bucket25_gate_audit.py` preserves its composer.
- `scripts/compose_koblitz_stage160_bucket27_gate_audit.py` and `scripts/compose_koblitz_stage161_filter8_gate_audit.py` preserve the later collector composers; their original audit directories already have collision-free `20260922` paths.

`GATE_STATUS.md` contains both narratives, with an explicit parallel-branch heading. `docs/index-calculus-scoreboard.html` retains the current-main scoreboard and adds the native-F4 single-target panel by element id. No scientific result or sealed native-F4 audit was recomputed to resolve the naming collision.
