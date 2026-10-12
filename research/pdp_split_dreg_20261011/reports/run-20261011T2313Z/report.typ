#set page(paper: "a4", margin: 16mm)
#set text(size: 9.5pt)

= Run report: PDP split/hybrid grid — amended grid (Amendment 1)

*Study* `research/pdp_split_dreg_20261011` (PROTOCOL.md). *Descriptive only*; interpretation is in the study README.

== What ran
- *Command:* `run_grid.sh` — per cell `pdp_split_dreg --cells CELL --guesses 0,2,4 --targets 4 --seed 20261011 --budget-secs 300 --max-free 24 --lambda-targets 256`, one process per cell under `timeout 7200`.
- *Records:* `runs-v2/*.jsonl` (one line per cell header and per run), stderr logs `runs-v2/*.stderr`, cell exits `runs-v2/exits.txt`.
- *Commit* `06d87db9f214f98233adae97778c772bf466b732`; *binary* SHA-256 `3752ee47cab471bc…` (`runs-v2/BINARY.sha256`).
- *When:* 2026-10-11 23:13:12 – 2026-10-12 00:20:05 UTC. *Host:* the same cloud container as the rho session (Intel Xeon @ 2.10 GHz, 4 vCPU); no isolation, so F4 wall time is descriptive only.
- *Cells* (m:a:n:ℓ): 3:1:11:4, 3:0:13:4, 3:1:17:4, 3:1:19:4, 3:0:13:5, 3:1:17:5, 4:1:11:3, 4:0:13:3. Plans unsplit and S₃ chain (m = 3), two-way 2+2 and three-way 2+1+1 (m = 4); `g ∈ {0, 2, 4}`; 4 targets per cell; 192 runs requested, 192 recorded.

== Outcome
All 192 runs complete and exact against exhaustive enumeration; no branch timed out; no false or unliftable solutions accepted. Against the first grid: solving degree and word XORs identical on all 192 runs; verified counts differ only on the three amended runs (0 → 1, 3, 1).

== Measured values
Per-run solving degree, F4 word XORs (deterministic), F4 wall (descriptive), verified and exhaustive decomposition counts are tabulated in `ANALYSIS.json` (and its Markdown table); the figures below plot those recorded values only.

#figure(image("fig1-solving-degree.svg", width: 100%), caption: [F4 solving degree per target, by plan and `g` (records `runs-v2/`, via `ANALYSIS.json`).])
#figure(image("fig2-word-xors.svg", width: 100%), caption: [F4 word XORs per target, all `2^g` branches summed, log scale (records `runs-v2/`, via `ANALYSIS.json`).])
#figure(image("fig3-projection-vs-rho.svg", width: 100%), caption: [Relation-phase projection in S at the pinned XOR ratio (derived, PROTOCOL.md formula) against measured rho from session `ECBS1he23809badf5b`; cells with a failed gate or `λ = 0` carry no projection.])
