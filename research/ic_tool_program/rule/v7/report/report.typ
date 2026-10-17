#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show raw: set text(size: 8pt)
= The rule's comparison at baseline v7

One unseen public point a process, the index calculus against Pollard rho on
that same point, at ledger §23's six sizes, 64 targets and four A/A repeats a
size: 408 rows, every one checked, every logarithm §23's. v7 is v6 with R11's
streamed pair-table build. Run on 2026-10-16 on an Emerald Rapids host of the
reference class. Protocol: `research/ic_tool_program/rule/v7/PROTOCOL.md`.

- *Verdict (§23's):* online the index calculus is faster than rho at all six
  sizes (14.1–37.0×). Cold, rho is still faster at all six.
- *Distance remaining:* rho is 3.43–11.51× faster cold (rule v6: 3.78–12.06×;
  rule v3: 5.23–20.97×). At the three largest sizes it is 7.03, 6.37 and
  11.51× (v6: 7.90, 7.40, 12.06×; v3: 13.24, 12.03, 20.97×). The gap grows as
  r^0.144 [0.092, 0.196].
- *Checks:* 408 of 408 rows pass the `vs_rho` claim check
  (verified_answer_controlled_wall); the pin held; the reference check is
  admissible at every size (rho at 0.57–0.78 of the strong walk's time a
  step); five rows' contended first attempts re-ran clean.
- *Pricing (AGENTS.md §5a):* both arms are timed in one process, so no native
  work is unpriced on either arm, and the 0 / 0.1 / 1 GAE band gives one
  figure.

#figure(image("rho-distance.svg", width: 100%), caption: [Left: cold time,
IC over rho (above 1, rho faster). Right: online speedup, rho over IC. Rules
v3, v6 and v7; 95% bootstrap intervals over 64 targets.])

#table(columns: 7, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*cold, v7 [95%]*], [*cold, v6*], [*cold, v3*], [*online, v7*], [*rho / strong*],
[`icv1-f2m47-t22705043-f4e44623`], [36.64], [3.43 [3.16, 3.75]], [3.78], [5.23], [22.59], [0.78],
[`icv1-f2m57-tm747311035-c1f545af`], [38.0], [3.97 [3.68, 4.33]], [4.54], [6.12], [20.87], [0.57],
[`icv1-f2m41-tm2308219-7f48b14a`], [39.0], [4.74 [4.27, 5.32]], [5.12], [7.66], [14.13], [0.76],
[`icv1-f2m53-tm56619371-dac20a85`], [44.26], [7.03 [6.27, 7.93]], [7.9], [13.24], [30.66], [0.72],
[`icv1-f2m59-tm943548413-98844ecc`], [44.52], [6.37 [5.6, 7.33]], [7.4], [12.03], [36.97], [0.73],
[`icv1-f2m61-t158598901-ab42b6c5`], [47.21], [11.51 [10.24, 13.06]], [12.06], [20.97], [19.46], [0.73],
)
