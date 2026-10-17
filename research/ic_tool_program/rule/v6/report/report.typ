#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= The rule's comparison at baseline v6

One unseen public point a process, the index calculus against Pollard rho on
that same point, at ledger §23's six sizes, 64 targets and four A/A repeats a
size: 408 rows, every one checked, every logarithm §23's. v6 is v4 (main's
`92246530`) with R06's scan key and R08's folding subtraction. Run on
2026-10-16 on an Emerald Rapids host of the reference class. Protocol and
amendment 1: `research/ic_tool_program/rule/v6/PROTOCOL.md`.

- *Verdict (§23's):* online the index calculus is faster than rho at all six
  sizes (14.2–36.4×). Cold, rho is still faster at all six.
- *Distance remaining:* rho is 3.78–12.06× faster cold (rule v3: 5.23–20.97×).
  At the three largest sizes it is 7.90, 7.40 and 12.06× (13.24, 12.03,
  20.97×). The gap still grows with the group, as r^0.142 [0.105, 0.178]
  (r^0.176 at v3).
- *Checks:* 408 of 408 rows pass the `vs_rho` claim check
  (verified_answer_controlled_wall). The pin held. The reference check is
  admissible at every size: rho at 0.70–0.85 of the strong walk's time a
  step. Four rows' contended first attempts re-ran clean.
- *Amendment 1:* the claim writer gained the ledger's five newer `vs_rho`
  fields, and the strong-rho fixture's re-imported size gate (which refused
  n = 47, 57 and 59) was widened. No row was re-timed.
- *Pricing (AGENTS.md §5a):* both arms are timed in one process, so no
  native work is unpriced on either arm, and the 0 / 0.1 / 1 GAE band gives
  one figure. The counted-unit readings stay lower bounds.

#figure(image("rho-distance.svg", width: 100%), caption: [Left: cold time,
IC over rho (above 1, rho faster). Right: online speedup, rho over IC. Rule
v3 grey, rule v6 blue; 95% bootstrap intervals over 64 targets.])

#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*cold, v6 [95%]*], [*cold, v3*], [*online, v6 [95%]*], [*online, v3*], [*rho / strong*],
[`icv1-f2m47-t22705043-f4e44623`], [36.64], [3.78 [3.48, 4.13]], [5.23], [21.97 [17.62, 28.07]], [10.17], [0.76],
[`icv1-f2m57-tm747311035-c1f545af`], [38.0], [4.54 [4.19, 4.95]], [6.12], [20.45 [16.3, 26.08]], [16.3], [0.85],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [5.12 [4.63, 5.72]], [7.66], [14.22 [11.3, 18.34]], [10.93], [0.75],
[`icv1-f2m53-tm56619371-dac20a85`], [44.26], [7.9 [7.03, 8.94]], [13.24], [30.66 [23.92, 39.61]], [21.47], [0.83],
[`icv1-f2m59-tm943548413-98844ecc`], [44.52], [7.4 [6.49, 8.52]], [12.03], [36.37 [28.18, 47.37]], [26.5], [0.7],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.21], [12.06 [10.72, 13.69]], [20.97], [19.88 [15.29, 26.51]], [14.82], [0.72],
)
