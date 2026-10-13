# Frontier synthesis replay report

Run id: `run-20261013T2025Z`. This was a native ECBench verification and
evidence-synthesis run, not a new timed benchmark. Operation-count replays are
the only new execution evidence; wall times in the frontier come from the
source studies and retain their original isolation labels.

## Build and host

- Repository commit: `759978250c1ce53d7b87bd8a1c2304781fcc4fc4`.
- Rust: 1.94.1, release profile, built in the official Bookworm container with
  the managed environment's CA bundle and TLS verification enabled.
- ECBench binary SHA-256:
  `53654ef1c0e93ed02025acd8188514d2674c4da2f57f25f61db40a510a422950`.
- Environment class: `ECBENV2h217559836c7c`; hardware class from the receipts:
  `ECBHW1h7373da65e433`.
- CPU: AMD EPYC 9V74, five logical CPUs exposed, one NUMA node, virtualized;
  Linux x86_64 6.18.44.

## Commands and outcomes

```text
ecbench verify factor-base K17 --replay 8
  30/30 records valid; 8/8 identical; receipt 9306d33e…

ecbench verify F4/F5/F6 K13 --replay 7
  48/48 records valid; 7/7 identical; receipt 0304c496…
  selected rows cover rho strong, rho negation, matrix F4, matrix F5,
  inherited F4, and F6-IC.

ecbench verify rootfinder K13 --replay 8
  224/224 records valid; 8/8 identical; receipt 4814136e…
  the deterministic sampler selected rho rows only.

ecbench verify rootfinder K13 --replay-all --exit-code
  FAILED: 120/168 identical, 48 inherited-F4 rows differ in solver phase
  counts; receipt 6110727a…. Every rootfinder and rho row is identical.
```

The failed receipt is retained in
[`../../audits/rootfind-n13-replay-all-failed.json`](../../audits/rootfind-n13-replay-all-failed.json).
The other receipts are beside it. No session directory was edited.

## Values rendered

- [`cold_ops_gap.svg`](cold_ops_gap.svg) plots committed cold `S / matched rho`
  ratios; hatching identifies lower bounds with unpriced work.
- [`online_boundary.svg`](online_boundary.svg) plots `rho / IC` on one common
  orientation. The n61 online and cold values are admitted; N83 online is
  exploratory, and its lower bar is only a setup-inclusive diagnostic—not a
  matched cold ratio. PNG copies are rendered from the same source at 150 dpi.
- [`../../frontier-data.json`](../../frontier-data.json) contains every plotted
  value, source path, and derivation.

## Observation limits

This run did not establish a new speedup, fit a new bound, or run a scaling
panel. It verifies representative current-binary coverage and exposes one
cross-host solver-counter replay failure. The result interpretation and exact
frontier definition are in [`../../README.md`](../../README.md).
