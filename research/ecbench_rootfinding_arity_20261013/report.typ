#set page(paper: "us-letter", margin: 0.72in)
#set text(font: "Helvetica Neue", size: 9.5pt)

= Prime pair-table arity at 64 abscissae

*Measured result.* Frozen spec commit `fd85dff96`, registered curve
`icv1-fp20-t727-cd198a38`, subgroup order 562427. Eight fresh public
targets, three measured rounds per arm, same point and algorithm seed in
each comparison. Session `ECBS1h4c7b9db8a389` is complete and sealed.

#table(
  columns: (1.7fr, 0.8fr, 0.9fr, 1.0fr),
  inset: 5pt,
  [Arm], [Verified], [Mean GAE], [Relation trials],
  [rho.negation], [24/24], [894.167], [—],
  [pair table m=2], [24/24], [17,909.292], [10,223],
  [pair table m=3], [24/24], [25,634.375], [227.5],
  [rootfind-subtract m=3], [0/24], [—], [setup error],
)

#image("cost.svg", width: 100%)

The paired m3/m2 counted quotient is *1.431345* (95% target-cluster
interval 1.224624–1.716117, 24 pairs). Mean setup is 4,160 GAE for
both pair-table arms; relation work is 12,087.917 GAE for m2 and 19,813
GAE for m3. The matched m2/rho counted quotient is 20.029; rho A/A is
1.000. Full local replay reproduced all 96 measured verified executions
with zero integrity problems.

The 32 errors belong to an auxiliary arm: the retained binary lacked the
new `rootfind-subtract` oracle and rejected it at setup. They remain
recorded, with no performance inference. Prime native work is partly
unpriced in GAE; every wall observation earned L0, so no runtime
comparison is admitted.

*Evidence:* `REPORT.md`, `summary.json`, `replay-all.json`,
`sessions/prime20-b64-m2-vs-m3/`, and `render_chart.mjs`. The figure is
this session's local counted-cost graph; no canonical frontier changed.
