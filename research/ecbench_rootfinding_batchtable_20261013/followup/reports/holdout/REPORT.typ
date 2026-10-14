#set page(margin: 18mm)
#set text(size: 10pt)

= Target-pinned factor-base holdout

_ecbench_ session `ECBS1hac7a406a30ba` on the registered 20-bit prime
subgroup used eight public targets unseen in the sweep. All 160 executions
verified and all 120 measured executions replayed identically. Audit
receipt SHA-256:
`3fb8c4b82a9e114e6537ec1aa459228c4ef08d142b3cf05cbe151425e801a7e4`.

#image("counts.svg", width: 100%)

B48 used 6,980.417 mean cold counted GAE versus 7,970.500 for B64:
paired ratio 0.8758, 95% target-cluster interval 0.8191–0.9370. B40
used 7,613.333; B48/B40 was 0.9169 with interval 0.8225–1.0150.
The B48 phase split was 2,352 table setup, 4,603.167 relations and
25.25 verification GAE. Matched negation rho used 1,047.000 mean GAE;
B48/rho was 6.6671 (5.3684–8.3883).

These are L0, lower-bound operation counts on one field size. Unpriced
native work prevents a full-cost claim. No wall comparison is admitted.

Source: `sessions/prime20-pinned-holdout-20261013/records.jsonl` and
`replay-all-holdout.json`. Figure source: `plot_counts.py`; plotted data:
`reports/holdout/summary.json`.
