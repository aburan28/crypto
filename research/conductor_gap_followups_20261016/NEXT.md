# Continue the authorized experiments

The user requested all five proposed experiments, a PR and tracking. Do not stop at preparation or silently substitute smaller cohorts. The serial campaign has 2,216 jobs; only immutable execution receipts establish actual completion. Live progress is in results/execution.jsonl and results/receipts. The parent launched as PID 17400. Inspect actual process state before acting; never launch a duplicate or remove its lock while alive.

## While running

- Preserve prepared/freeze.json and its listed files/executables. The benchmark source revision is in prepared/source-commit.txt. New post-hoc source files are separate from the frozen benchmark.
- Do not compile or run other solvers while a measured child is active. For necessary maintenance, stop only the scheduler, let the complete Rust child (including replay/archive) and GNU time wrapper finish, verify no measured descendant remains, then work and resume. A native-backend exit alone is insufficient; see results/timing-deviations.json.
- Use target/release/examples/conductor_gap_followups progress for lightweight metadata. No full bootstrap while the campaign lock exists. Queue order is crossover, yield, rho, coordinate controls, orbit audit, dose response and overlap-selected/held-out solves.
- Active files are deliberately ignored. Publish checkpoints by staging immutable receipt/result directories and their referenced compressed base artifacts; do not stage a current writer or claim a queued panel is complete.
- Keep all failures, timeouts and censors. The scheduler preserves receipts and never overwrites a run. If a panel fails, name the actual gap; retain the full original workload rather than successful-only filtering.

## After the final job

1. Confirm all 2,216 receipts and final campaign event. The lock must be absent. Preserve the entire final execution journal and host/resource snapshots.
2. Run /workspace/native-deps/cm-followups-v1 audit-all. It independently checks points, rows/rank, query indices and recovery, including all censored observations. Execution failures fail closed and need an explicitly declared repair/rerun panel.
3. Run target/release/examples/conductor_gap_followup_posthoc analyze. This separate postprocessor checks common query prefixes for wall-censored same-crater arms; do not use the frozen driver's stricter analyze command. Its regressions reject changes inside the common prefix and unexplained unequal query counts. Inspect counts, all censors and the recorded timing deviation before interpreting ratios. Add a CPU sensitivity analysis excluding the entire gap affected by crossover-0003; retain the full primary result as well.
4. Complete the five per-requirement result reports, rank/exposure summaries, dose-response class contrasts and overlap-selection interpretation. The orbit panel is arithmetic-only; recoveries may pin the target at rank four without solving every seed log. A held-out comparison is not an independent population confirmation. The dose panel has four classes per CM field and small conductor products, not 96-bit ECDLP benchmarks.
5. Generate independent PARI target replay vectors from sealed successful records if practical, and preserve all verification receipts. Export a PDF summary and SVG/150-dpi PNG figures from native-derived metrics. Keep rendering separate from arithmetic/statistical analysis; no Python research harness.
6. Register the 32 newly certified models/identities and update the affected canonical catalogs, study/dashboard/browser panels and benchmark/problem indexes under AGENTS.md. Do not invent eligible speedup scores for unpriced native units. No complete IC/rho cost ratio is established by these generic MITM comparisons.
7. Run relevant validation, update evidence manifests, commit and push to existing PR #1888. Preserve the original manifest snapshot under inputs/. Do not modify frozen historical source pins or archives to hide inherited CI failures. Required human review and applicable checks still gate merge.

## Tracking

A one-time automatic check-in was created to resume this task. Re-arm under the repository cadence while work advances; do not send repeated no-change updates. Check actual access rather than claiming a missing environment/process is running. Only three neutral/success checks and one LLM workflow were visible on head 80298bc8; applicable native/standard CI was not established green and PR mergeability was unknown. Recheck on a substantive new head, rather than pushing empty commits. The CLI and managed runtime had working read access; conductor/cairn binaries were absent and no distributed lease/funding action is claimed.
