#set page(paper: "us-letter", margin: 0.72in)
#set text(font: "Helvetica Neue", size: 9.5pt)

= Ecbench run: prime pair-table arity

Session `ECBS1h4c7b9db8a389`, started 2026-10-13 17:06 UTC. Frozen spec
`specs/prime20-b64-m2-vs-m3.json` was committed at `fd85dff96` before
execution. The binary SHA-256 is `be956553a5a2cd1227e8f1add8dd24e11dffba5068c799234b55e44ffc3a1b61`.
The host is macOS arm64; `--cpus none --allow-busy` gave L0 isolation.
Eight public targets of `icv1-fp20-t727-cd198a38` were run with target
seed 20261014, one warmup and three measured rounds per arm.

#table(
  columns: (1.7fr, 0.8fr, 0.9fr, 1.0fr),
  inset: 5pt,
  [Arm], [Verified], [Mean GAE], [Mean relation trials],
  [rho.negation], [24/24], [894.167], [—],
  [pair table m=2], [24/24], [17,909.292], [10,223],
  [pair table m=3], [24/24], [25,634.375], [227.5],
  [rootfind-subtract m=3], [0/24], [—], [setup error],
  [rho A/A], [24/24], [894.167], [—],
)

#image("../../cost.svg", width: 100%)

The auxiliary rootfind arm emitted 32 setup errors because this retained
binary does not register that oracle on prime curves. The session retains
all errors. The audit `../../replay-all.json` accepted 160 records and
reproduced 96 measured verified runs with zero problems. Raw records are
in `../../sessions/prime20-b64-m2-vs-m3/records.jsonl`; exact run and
verify commands are in `../../REPORT.md`. Prime native costs are not
fully priced in GAE; the L0 wall values are descriptive only.
