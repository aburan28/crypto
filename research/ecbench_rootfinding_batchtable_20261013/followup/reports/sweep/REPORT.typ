#set page(margin: 18mm)
#set text(size: 10pt)

= Target-pinned relation collection: frozen sweep

_ecbench_ session `ECBS1hdabf729ffad4` on registered prime curve
`icv1-fp20-t727-cd198a38` (order 562427). Eight public targets, three
measured rounds, one warmup. All 256 executions verified; all 192 measured
executions replayed identically. The exact spec, session and audit receipt
are stored in this study. Audit receipt SHA-256:
`1a24f267ecec04dd337c2c64f603273f8a48e4058ab4c2a67d84d991c2e04dd1`.

#image("counts.svg", width: 100%)

The bars include per-target table setup, relation collection and exact
verification in counted group-addition equivalents (GAE). B48 target pin
had the lowest mean in the prespecified B40/B48/B64/B80 sweep at 7,290.083
GAE. B64 target pin used 8,082.625 GAE versus 18,740.917 for checked full
rank (paired ratio 0.4313, 95% target-cluster interval 0.3987–0.4701).
B48 / B64 was 0.9019 (0.8591–0.9488); B48 / B40 was 0.9595
(0.8853–1.0343). The current-binary B64 subtraction arm verified but
used 342,765.167 mean GAE, 42.4077 times the B64 table arm.

All executions were L0 on shared macOS hardware. Unpriced native work
means the counted ratios do not bound full cost. The best IC arm remained
7.1962 times matched negation rho in counted GAE at this curve size.

Source: `sessions/prime20-pinned-sweep-20261013/records.jsonl` and
`replay-all-sweep.json`. Figure source: `plot_counts.py`; plotted data:
`reports/sweep/summary.json`.
