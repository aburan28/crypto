#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
= LLL degeneracy probes, re-run 2026-10-16

This report records observations only. It re-runs the committed probes behind
`papers/lll-degeneracy-note` to correct that note (sheriff case
SHF-20261016-3caa6f, blocker P1). Code: the commit in `logs/commit.txt`.
Commands and logs are in `README.md` and `logs/`.

#figure(image("recovery.svg", width: 100%), caption: [Share of probes recovering the planted key: the note's pre-fix reports (2026-10-11) against the 2026-10-16 re-run (`results.json`).])

- *secp256k1:* 6/6 recovered. The note reported 0/3 and called it "curve-specific".
- *P-384 and brainpoolP384r1:* 4/4 each. The note reported iteration-cap failures.
- *P-521:* 0/4 on the scaled `f64` path, and 3/3 with 2048-bit fixed-point Gram–Schmidt.
- *Every 256- and 384-bit probe* finished in 13–27 ms (one unisolated 4-vCPU VM; not a benchmark).

Both failure modes were the `f64` Gram–Schmidt bugs that
`research/notes/lattice-hnp/RESEARCH_LLL_GS_ANALYSIS.md` §10 recorded: an
overflow of $n^4$ (secp256k1, all 384-bit curves) and catastrophic cancellation
(P-521).

*Limits:* 1 to 6 probes per curve at $m = 8$. These results show the reported
failures do not reproduce. They do not make a security claim or a cost claim.
