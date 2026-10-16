#set page(paper: "a4", margin: 2cm)
#set text(size: 9.5pt)
= IC frontier campaign — Phase 1 run report

*Run.* 2026-10-15, sessions 07:09–08:43 UTC. Protocol `PROTOCOL.md` (frozen before the run), specs
`phase1/spec-*.json`, parameters `params.json`. Host: repository cloud container, Linux x86_64, 4 vCPU,
15 GiB; `--cpus none --allow-busy`: every record is *L0*, operation counts only, no wall time admitted.
The `ecbench` binary was rebuilt during session P32 with analysis-only changes (`phase1/BINARY_NOTE.txt`).

*Workloads.* Public hash-to-subgroup hold-out targets (seed 20261015200), 16 per curve, 2 measured rounds
and 1 warm-up. Koblitz $(a,n)$ ∈ {(0,13), (1,17), (1,19), (1,23), (0,37), (1,43), (1,47), (0,41)};
$a=-3$ prime-order curves of 20, 24, 28, 32, 36, 40 bits; PKM ladder curves of 20–28 bits.

*Requested and completed.* 19 sessions, 4,848 executions, 3,232 measured: 3,178 verified and 54 exhausted
(all `ic-shared`: 30 at $n=41$, 24 at $n=47$). No timeouts, errors or crashes. Replay-all audits:
19/19 sessions `ok`, every replayed record identical (`phase1/audits/`).

#figure(image("frontier.svg", width: 100%),
  caption: [Mean $S$ per arm and curve at $w = 0.1$ GAE per native operation, bar from $w = 0$ to $w = 1$
  (`ecbench band`, `analysis/band-*.json`). Open marker: not every run verified.])

*Observations (the harness's band verdicts, `analysis/band-*.md`).* On every Koblitz curve
`rho-frob` < `ic-compact` and `rho-frob` < `ic-negfold` (where it runs, $n <= 37$) at all three weights; the one unresolved
IC–rho pair is `ic-compact` vs `rho-strong` at $n = 13$. On every prime and PKM curve `rho-neg` < every
IC arm at all three weights. IC/best-rho mean-S ratios: `analysis/ratios_to_rho.tsv`. Exponent fits per
arm, tier and weight (`ecbench bound fit --weight [--param-rule]`): `analysis/bounds/` (63 records);
no IC arm has a fitted exponent below the matching rho exponent. `ic-shared` has uncounted unpriced categories
and is reported at $w = 0$ only.

*Raw records.* `phase1/sessions/<curve>/records.jsonl`. Interpretation: `README.md` and the paper.
