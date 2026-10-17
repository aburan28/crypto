#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
= LLL degeneracy probes, re-run 2026-10-16

This report records observations only. It re-runs the committed probes behind
`papers/lll-degeneracy-note` to correct that note (sheriff case
SHF-20261016-3caa6f, blocker P1). Code: the commit in `logs/commit.txt`.
Commands and logs are in `README.md` and `logs/`.

#block(stroke: 0.6pt, inset: 6pt, width: 100%)[*Correction (2026-10-17, sheriff case SHF-20261017-e4d7f6, P-N1 and P-N2).* Counts were probe runs of deterministic, overlapping tests; distinct seed pairs are now given. Every P-521 probe ran at $k_"bits" = 384$, not 390 (= ⌊3·521/4⌋). The attribution of the P-521 scaled-`f64` failure to catastrophic cancellation is withdrawn as unsupported (with the scaling, P-521's diagonal is exactly $2^500$ in `f64`); the measured results stand. The figure counts runs; every share is 0 or 1 either way.]

#figure(image("recovery.svg", width: 100%), caption: [Share of probes recovering the planted key: the note's pre-fix reports (2026-10-11) against the 2026-10-16 re-run (`results.json`).])

- *secp256k1:* 3/3 distinct seed pairs recovered (6/6 probe runs). The note reported 0/3 and called it "curve-specific".
- *P-384 and brainpoolP384r1:* 3/3 distinct seed pairs each (4/4 runs). The note reported iteration-cap failures.
- *P-521* ($k_"bits" = 384$, not 390): 0/3 distinct seed pairs (0/4 runs) on the scaled `f64` path, and 3/3 with 2048-bit fixed-point Gram–Schmidt.
- *Every 256- and 384-bit probe* finished in 13–27 ms (one unisolated 4-vCPU VM; not a benchmark).

Both failure modes were the `f64` Gram–Schmidt bugs that
`research/notes/lattice-hnp/RESEARCH_LLL_GS_ANALYSIS.md` §10 recorded: an
overflow of $n^4$ (secp256k1, all 384-bit curves) #strike[and catastrophic cancellation
(P-521)]. P-521 fails on the scaled `f64` path for a reason not yet isolated.

*Limits:* 1 to 3 distinct seed pairs (1 to 6 probe runs) per curve at $m = 8$. These results show the reported
failures do not reproduce. They do not make a security claim or a cost claim.
