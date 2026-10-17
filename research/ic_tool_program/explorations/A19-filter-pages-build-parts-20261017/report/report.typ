#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
= A19: the huge-page filter with the folded build's 64 partitions

*Exploration, no round.* v8 (`ic-r12-v7-23138d20`) against A17's huge-page filter plus the folded build's 64 partitions (`ic-a19-v8-6c998fb0`, `candidate.patch`). Run 2026-10-17 06:44–07:02 UTC: 72 isolated processes, none contended; every logarithm the same in both arms.

*Outcome:*
- *2^44.3:* 1.067.
- *2^47.2:* 1.054.
- *2^44.5:* 1.014 [0.978, 1.051]. Its collection does not respond, as in A16 and A17.
- *No round.* A round accepts only if each size's interval has its lower end above 1.03 (`icprog`'s `bound_reasons`), and 2^44.5's lower end is below 1.

#figure(image("cold-ratio.svg", width: 100%), caption: [Cold time, v8 over a19, per pair, with the geometric mean and 95% t interval. Dashed: 1.03.])

#table(columns: 6, stroke: 0.4pt, inset: 3pt,
  [*curve*], [*log₂ r*], [*v8 / a19 [95%]*], [*collection*], [*build*], [*selection*],
  [`icv1-f2m53-tm56619371-dac20a85`], [44.3], [1.067 [1.010, 1.126]], [1.045], [1.096], [1.091],
  [`icv1-f2m59-tm943548413-98844ecc`], [44.5], [1.014 [0.978, 1.051]], [1.003], [1.048], [1.025],
  [`icv1-f2m61-t158598901-ab42b6c5`], [47.2], [1.054 [1.012, 1.097]], [1.071], [1.014], [1.015],
)

*The build-timing diagnostic* (`partbits-sweep.tsv`; single processes; not a measurement): per stored pair, add and key cost about 8 ns, the append into the streams 12–13 ns, and the place step 9.5–11.6 ns. With 64 partitions the append is about 30% cheaper.

Scope: one host, one core, three sizes, two rows a size.
