#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9pt)
= Run report: run-20261014T1925Z (ecbench, Frobenius-orbit IC vs rho, fully priced)

*Command.* `bash research/ecbench_frobenius_orbit_vs_rho_20261014/run.sh` (eight `ecbench run --wait --cpus none --allow-busy`, then `verify --replay-all`, `table`, `compare --save`). Specs `docs/ecbench/specs/claim-koblitz-frob-orbit-full-a<a>-n<n>.json`, frozen in `PROTOCOL.md` at commit `e1708eef8`; binary sha256 `f9cfb1a09cc65074b03ee4af1aab71438ab11b2467391abcfed1594a98484d19` built from that commit. Start 2026-10-14T19:25:54Z (UTC).

*Host and isolation.* Intel Xeon @ 2.80 GHz, 4 vCPU, shared container; L0 on every record (no CPU reservation). Wall times are descriptive only.

*Workloads.* Eight Koblitz curves, 16 public targets each (target seeds 202610142001–202610142015, odd), 3 measured rounds + 1 warm-up, alternate order. Arms: rho-frob (`rho.signed_frobenius_full`, reference), rho-frob-aa (A/A control), rho-strong (`rho.signed_frobenius_strong_full`), `ic-orbit-c<c>-renew` and `-walk` (`ic.pipeline_full`, compact-orbit-scan, mitm-frobenius-counted:m=2), `ic-negfold-d<d>-renew` (n ≤ 23 only).

*Completion.* 2,816 records requested, 2,816 written, 2,816 verified; no failures, no timeouts. `verify --replay-all`: 2,112 of 2,112 measured records replayed identically; every audit `ok`. Two sessions (a1-n17, a1-n29) record `git_dirty = true` (whitespace-only rustfmt of source during the run; binary unchanged).

*Measured values* (mean S = GAE/√r over 48 measured runs per arm, 95 % interval as `ecbench table` recorded; every record `lower_bound = false`). Unit: `ecbench.gae`.

#set table(inset: 4pt)
#table(columns: (10.5em, 3.6em, 1fr, 1fr, 1fr, 1fr, 1fr), [*session*], [*curve*], [*rho-frob*], [*rho-strong*], [*orbit renew*], [*orbit walk*], [*negfold renew*],
[#text(size: 7.5pt)[`ECBS1hb67a421a6fe2`]], [a0-n13], [2.025 [1.863, 2.321]], [14.562 [14.255, 15.026]], [5.986 [5.857, 6.116]], [14.847 [14.672, 15.052]], [11.718 [10.688, 12.801]],
[#text(size: 7.5pt)[`ECBS1h2627f1e2aee2`]], [a1-n17], [1.524 [0.921, 3.218]], [4.153 [4.026, 4.310]], [2.619 [2.348, 2.932]], [5.210 [4.867, 5.586]], [11.214 [10.501, 11.972]],
[#text(size: 7.5pt)[`ECBS1h8baa27b989e6`]], [a0-n19], [0.835 [0.733, 1.009]], [3.172 [3.084, 3.262]], [2.343 [2.233, 2.478]], [4.217 [4.100, 4.362]], [15.786 [13.819, 17.902]],
[#text(size: 7.5pt)[`ECBS1h5ace8f3e67d1`]], [a0-n23], [0.416 [0.344, 0.516]], [1.161 [1.109, 1.215]], [1.601 [1.468, 1.742]], [2.126 [1.992, 2.266]], [24.807 [23.361, 26.152]],
[#text(size: 7.5pt)[`ECBS1h575c8bcdc35a`]], [a1-n29], [1.665 [1.109, 2.639]], [4.543 [4.454, 4.644]], [2.607 [2.532, 2.699]], [5.431 [5.346, 5.528]], [–],
[#text(size: 7.5pt)[`ECBS1ha274653650f8`]], [a0-n31], [0.394 [0.330, 0.485]], [1.295 [1.245, 1.350]], [1.445 [1.342, 1.555]], [2.163 [2.055, 2.286]], [–],
[#text(size: 7.5pt)[`ECBS1h7aed5aa2d512`]], [a0-n37], [0.253 [0.217, 0.290]], [0.273 [0.248, 0.298]], [1.275 [1.153, 1.410]], [1.239 [1.127, 1.350]], [–],
[#text(size: 7.5pt)[`ECBS1h868a2da9acdc`]], [a0-n41], [0.133 [0.113, 0.154]], [0.157 [0.132, 0.184]], [2.412 [2.301, 2.526]], [2.390 [2.267, 2.523]], [–],
)

#figure(image("s-by-arm.svg", width: 100%), caption: [Mean S per arm and curve with the recorded 95 % intervals; sessions as in the table.])

*Controls.* rho-frob-aa ÷ rho-frob = 1.000 [1.000, 1.000] on every curve (`tables/compare-*-rho-frob-aa-vs-rho-frob.txt`).

*Raw records.* `sessions/<curve>/records.jsonl`; audits `audits/audit-<curve>.json`; tables and comparisons `tables/`, `sessions/<curve>/comparisons/`. Interpretation is in the study README and REPORT.pdf.
