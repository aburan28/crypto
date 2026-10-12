#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261011T2258Z-degree-bound]

*What ran.* `cargo run --release --example isogeny_loop_degree_bound -- research/isogeny_loop_glv_20261010/group_{a_koblitz_and_known,b_binary131,c_certicom_mid,d_certicom_small}.json --json research/isogeny_loop_glv_20261011/degree_bound.json`. Base commit db0393b4 plus the example source with SHA-256 `3b93c9f1…d27b` (see `../../revision.txt`); rustc 1.97.0; Linux 6.18 x86_64 cloud container, 4 vCPU, not isolated (no timing claim is made). Start 2026-10-11T22:58Z (UTC); wall time under one second. An earlier run at 22:47Z, from the same source before formatting and clippy-only edits, produced identical rows.

*Inputs.* The frozen round-1 sweep records (Python provenance, `research/isogeny_loop_glv_20261010/`): per curve $q$, trace $t$, subgroup order $r$, fundamental discriminant $D_K$, conductor $f$, cheapest loop found.

*Runs requested / completed.* 1 / 1, exit 0. 21 curves read; for all 21 the identity $t^2 - 4q = f^2 D_K$ re-verified in native big-integer arithmetic (column "verified": yes for every row); for every prime-field curve the bound was asserted not to exceed the cheapest loop found.

*Recorded values (derived per curve from the inputs; unit: M = field multiplications).* Bound $= kappa log_2 d_min (D_K)$ with $kappa$ = 14.196 / 5.047 / 3.786 (projective / affine / optimistic x-only); saving $= 8 ceil(log_2 r slash 2)$. Full table: `../../degree_bound.md`; machine record: `../../degree_bound.json`. Range over the 13 random challenge curves: bound/saving 3.27–3.51 (projective), 1.16–1.25 (affine), 0.87–0.94 (x-only); CryptoPro-A 3.55 / 1.26 / 0.95; CryptoPro-B 0.10 / 0.04 / 0.03; Koblitz curves and secp256k1 trivial ($d_min <= 2$).

#figure(image("../../figures/degree_bound.svg", width: 80%), caption: [Recorded bound ratios and round 1's cheapest loops per curve (run run-20261011T2258Z-degree-bound).])

*Observations only.* Interpretation is in `../../README.md` §1 and `../../REPORT.pdf`.
