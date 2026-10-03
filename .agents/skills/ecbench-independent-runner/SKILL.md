---
name: ecbench-independent-runner
description: Run an ecbench ECDLP measurement on an independent lab runner (isolab) with whole-core NUMA-local placement, and reproduce or cross-check a committed ecbench session on another machine with exact replays and a replay certificate. Use when a result needs L2/L3 isolation, a second host, or independent validation before a claim.
---

# Independent runs and reproduction

Read [`docs/ecbench/README.md`](../../../docs/ecbench/README.md) §6 and §9, and
`isolab/README.md` for the lab itself. ecbench writes the job and isolab places
and isolates it. Every run is still graded by ecbench from its own observations.

## Submit a spec to isolab

1. Build `ecbench` for the worker's architecture (Linux x86-64 or arm64), from a
   committed, clean tree. The binary's SHA-256 is its identity.
2. Write the job:

   ```bash
   ./target/release/ecbench isolab-job --spec SPEC.json --binary path/to/linux/ecbench --cpus 2 --policy strict > job.json
   ```

   - `--cpus 2` places one core and its SMT sibling. The measured child runs on
     one CPU and the runner sleeps on the other. Use `--cpus 4` for two cores.
   - `--policy strict` (the default) expands in isolab to bare metal, perf
     counters and isolation tier A, so it places only on a bare-metal worker. On
     a VM worker use `--policy standard` and expect the records to show steal or
     frequency blockers.
   - `--git-url URL --commit SHA` builds on the worker instead. `Cargo.lock` is
     untracked, so dependencies resolve there and the job needs network. Prefer
     `--binary`.
3. Submit with the isolab MCP tools (`isolab_submit`, then `isolab_wait`) or the
   CLI (`isolab submit job.json`). Do not edit the generated `fidelity` or
   `resources` sections to make a job place.
4. Fetch the artifacts (`isolab fetch JOB_ID`, or `isolab_result` over MCP).
   `session/` is an ordinary ecbench session. The job's
   `verify` step already ran `ecbench verify --replay 3` on the worker. Run the
   audit again locally:

   ```bash
   ./target/release/ecbench verify --dir artifacts/session --replay 12 --out artifacts/session.audit.json
   ```

5. Report the levels earned and the top blockers from the records, not from the
   job's policy. A strict placement that still shows steal or foreign ticks earned
   less, and the record says so.

If `meow34.pem` or lab credentials are not available in your environment, report
that blocker. Do not create replacement keys or accounts (AGENTS.md §9).

## Reproduce a session on another machine

Operation counts depend only on code, workload and seed, never on the host. Any
machine can therefore check them:

```bash
./target/release/ecbench verify --dir research/TOPIC/sessions/SESSION --replay 24 --exit-code --out receipt.json
```

- Use a binary built from the session's `git_commit` (in `session.json`). A
  different commit is a different method unless the replays still match.
- `identical` on every replay means the figure is reproduced. Cite the receipt's
  SHA-256 as the independent replay certificate.
- Wall time does not reproduce across hosts, and nothing in this check claims it
  does. To confirm a wall-clock result, rerun the spec as a new session on a host
  of the same class (`env_class_id`) at the required level and compare the two
  intervals.

## When results disagree

Do not average them away. Record both sessions. Compare `host.json` stable facts
(the class id says whether the hosts match), the levels and their blockers, and
`binary_sha256`. Report the disagreement as a finding.
