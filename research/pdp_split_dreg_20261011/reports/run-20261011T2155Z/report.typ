#set page(paper: "a4", margin: 16mm)
#set text(size: 9.5pt)

= Run report: PDP split/hybrid grid — first grid (pre-amendment)

*Study* `research/pdp_split_dreg_20261011` (PROTOCOL.md). *Descriptive only*; interpretation is in the study README.

== What ran
- *Command:* `run_grid.sh` — per cell `pdp_split_dreg --cells CELL --guesses 0,2,4 --targets 4 --seed 20261011 --budget-secs 300 --max-free 24 --lambda-targets 256`, one process per cell under `timeout 7200`.
- *Records:* `runs/*.jsonl` (one line per cell header and per run), stderr logs `runs/*.stderr`, cell exits `runs/exits.txt`.
- *Commit* `51bb1dc5e1fea5a7b53a2c6c8b02f51f75e77233`; *binary* SHA-256 `e8b4559915ef50a9…` (`runs/BINARY.sha256`).
- *When:* 2026-10-11 21:55:57 – 22:59:18 UTC. *Host:* the same cloud container as the rho session (Intel Xeon @ 2.10 GHz, 4 vCPU); no isolation, so F4 wall time is descriptive only.
- *Cells* (m:a:n:ℓ): 3:1:11:4, 3:0:13:4, 3:1:17:4, 3:1:19:4, 3:0:13:5, 3:1:17:5, 4:1:11:3, 4:0:13:3. Plans unsplit and S₃ chain (m = 3), two-way 2+2 and three-way 2+1+1 (m = 4); `g ∈ {0, 2, 4}`; 4 targets per cell; 192 runs requested, 192 recorded.

== Outcome
One cell failed correctness gate 2: `3:0:13:5`, plan `split-2+1`, `g = 0`, targets 0, 1 and 3 read out no decompositions (exhaustive: 1, 3, 1). Cause: solution read-out on a non-inter-reduced F4 basis (AMENDMENT_1.md). All other 189 runs complete and exact. A first launch of cell 3:1:11:4 was stopped by the operator to relaunch the orchestration detached; its 24 completed runs are kept in `runs/aborted-launch-1/` and reproduce this grid's solving degrees and XOR counts exactly.

== Measured values
Per-run solving degree, F4 word XORs (deterministic), F4 wall (descriptive), verified and exhaustive decomposition counts are tabulated in `ANALYSIS-runs-v1.json` (and its Markdown table); the figures below plot those recorded values only.

#figure(image("fig1-solving-degree.svg", width: 100%), caption: [F4 solving degree per target, by plan and `g` (records `runs/`, via `ANALYSIS-runs-v1.json`).])
#figure(image("fig2-word-xors.svg", width: 100%), caption: [F4 word XORs per target, all `2^g` branches summed, log scale (records `runs/`, via `ANALYSIS-runs-v1.json`).])
#figure(image("fig3-projection-vs-rho.svg", width: 100%), caption: [Relation-phase projection in S at the pinned XOR ratio (derived, PROTOCOL.md formula) against measured rho from session `ECBS1he23809badf5b`; cells with a failed gate or `λ = 0` carry no projection.])
