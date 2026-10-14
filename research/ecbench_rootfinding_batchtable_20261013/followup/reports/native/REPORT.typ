#set page(margin: 18mm)
#set text(size: 10pt)

= Native one-lane relation walk at B48

_ecbench_ session `ECBS1h8d99a8f4f164` used eight fresh public targets on
the registered 20-bit prime subgroup. All 128 executions verified, and all
96 measured executions replayed identically. The 24 paired native/generic
runs matched relation stream, matrix counts, target scalar and counted cost.
Binary SHA-256: `1ff5b6d5f717f79cdfc557fde0017b0fac7e40502fe26a28bd7cdb2a99354904`.

#image("counts.svg", width: 100%)

Both B48 arms averaged 7,869.292 cold counted GAE, with 2,352 table setup,
5,490.542 relations and 26.750 verification; ecbench's paired operation
ratio was 1.0000. The one-lane native path made 3,965.375 walk additions
per measured execution on average. The L0 median native/generic solve-wall
ratio was 1.3898, with no L2-admitted pair. These timings are descriptive.
The native arm is therefore a correctness baseline for subsequent lane
batching, not a demonstrated speedup.

Matched negation rho used 989.833 mean counted GAE; native/rho was 7.9501
(95% target-cluster interval 5.9335–10.8337), with native and rho costs
still incompletely priced.

The separate lane kernel preserved exact final points and reduced affine
inversions from up to one per lane to one per synchronized step at widths
4, 8 and 16. It excludes target probes and full-pipeline costs.

Sources: `sessions/prime20-native-walk-20261013/records.jsonl`,
`replay-all-native.json`, `native-pair-check.json`, and `native-lane-profile.txt`.
Figure source: `plot_counts.py`; plotted data: `reports/native/summary.json`.
