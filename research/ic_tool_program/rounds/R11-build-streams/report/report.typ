#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= R11: the pair-table build by partition streams and per-run filters

*Accepted, as baseline v7* (engineering). v6 (v4 plus R06's scan key and
R08's folding subtraction) against v6 plus R11's candidate
(`candidate-v6.patch`: each pair's `(hash, word)` straight into one of 256
partitions on its run's stream, and a presence filter per run ORed at the
end), both built with `rustc 1.97.0`, on an Emerald Rapids host of the
reference class, 2026-10-16. Protocol:
`research/ic_tool_program/rounds/R11-build-streams/PROTOCOL.md`
(declaration, amendment 1); analysis: `analysis.json` (`icprog analyse
r11`); run trees: `runs.tar.xz` (SHA-256 `5fe8adebbcc9e02f…`).

- *Tests:* 207 of 207 `koblitz_fast` and `koblitz_index_calculus` tests
  on the candidate's tree, among them the word-for-word table test at one,
  two, three and five threads.
- *Pin:* the candidate equals v0's outputs on all 90 suite rows.
- *Timed:* 946 isolated processes (A/A 224, suite 401, holdouts 321),
  6 contended and re-run, none failed; the driver task was relaunched once
  on the same container and resumed from disk.
- *Decision:* at each target size the cold interval lies above 1.03 on the
  suite rows and on the fresh holdouts; no size regresses beyond its A/A
  band. Predicted 1.08–1.12, 1.08–1.12 and 1.04–1.07.
- *The build phase* (medians over every repetition, base → candidate):
  187.1 → 124.5 ms, 190.4 → 126.3 ms and 331.1 → 219.0 ms, 1.50–1.51×
  faster. Collection is unchanged (0.98–1.01).

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v6 over the
candidate (above 1, the candidate is faster), paired 95% t intervals; grey:
the round's own A/A band; dashed: 1.03. Values from `analysis.json`.])

== Suite rows (five rounds)
#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*rows*], [*cold [95%]*], [*collection*], [*S cold, v6 → v7*], [*A/A band*],
[`icv1-f2m19-tm797-9c54981b`], [18.002], [2], [1.011 [0.972, 1.052]], [0.795], [124.09 → 122.43], [0.919–1.034],
[`icv1-f2m23-tm5197-1f85e9e1`], [22.001], [2], [1.021 [0.995, 1.047]], [1.005], [32.33 → 31.42], [0.938–1.154],
[`icv1-f2m45-tm6236725-40939294`], [24.803], [2], [1.019 [0.992, 1.047]], [0.982], [19.77 → 19.38], [0.879–1.187],
[`icv1-f2m37-tm534059-32aad96b`], [27.781], [2], [1.043 [0.962, 1.131]], [1.015], [8.91 → 8.44], [0.973–1.057],
[`icv1-f2m43-tm998717-e2e742b0`], [32.113], [2], [1.017 [0.967, 1.068]], [0.999], [4.57 → 4.5], [0.968–1.01],
[`icv1-f2m47-t22705043-f4e44623`], [36.636], [2], [1.061 [1.026, 1.098]], [1.003], [2.53 → 2.4], [0.988–1.054],
[`icv1-f2m57-tm747311035-c1f545af`], [38.002], [2], [1.137 [1.046, 1.236]], [1.049], [2.91 → 2.62], [0.909–1.02],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [2], [1.102 [1.041, 1.166]], [0.981], [3.38 → 3.1], [0.947–1.026],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [8], [1.134 [1.107, 1.162]], [1.014], [3.87 → 3.48], [0.936–1.026],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [8], [1.124 [1.1, 1.148]], [0.982], [3.53 → 3.11], [0.939–1.064],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [8], [1.071 [1.049, 1.093]], [1.011], [4.82 → 4.51], [0.986–1.056],
)

== Fresh holdouts (seeds 234–237, T159–T166; ten rounds at 2^44.3)
#table(columns: 5, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*pairs*], [*cold [95%]*], [*S cold, v6 → v7*],
[`icv1-f2m53-tm56619371-dac20a85`], [44.259], [80], [1.126 [1.104, 1.148]], [3.83 → 3.38],
[`icv1-f2m59-tm943548413-98844ecc`], [44.517], [40], [1.157 [1.129, 1.186]], [3.72 → 3.15],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.211], [40], [1.065 [1.048, 1.083]], [4.57 → 4.26],
)

== Scope
One host and one core, the reference class. The table is word for word the
same, by construction and by test, so the counts and the ratio to the floor
are unchanged. Chained with R06 and R08, cold time is 1.80–1.86× faster than
v4 at the three target sizes (a product of paired ratios, not a measurement).
The distance to matched rho is the rule at v6's until a rule runs at v7.
