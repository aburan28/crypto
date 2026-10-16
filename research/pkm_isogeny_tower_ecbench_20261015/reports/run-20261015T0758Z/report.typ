#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#show heading: set text(size: 11pt)

= Run report: pkm_isogeny_tower_ecbench_20261015, sessions p20–p28 and q20–q24

Descriptive record of the runs (observations only; interpretation is in the study README).

- *Started* 2026-10-15T07:58:18Z, finished 08:00:32Z (`progress.txt`); panel 3 (`run_panel3.sh`) after.
- *Binary*: `ecbench` built from commit `43986439e`, SHA-256 `074288d4a5b88b683f851d0a49900ea7b016a0f4c66651fb01cda16c748a13d3`.
- *Host*: Intel Xeon @ 2.10 GHz, 4 vCPU, virtual, shared with four other jobs; env class `ECBENV2ha95bed7071ef`. `ecbench run --cpus none --allow-busy`: every run L0. One thread, `nice -n 10`, 2.5 GB address-space cap.
- *Specs*: `specs/SPEC-{p20,p22,p24,p26,p28,q20,q22,q24}.json`; 8 one-target workloads per curve, 2 measured rounds + 1 warm-up, order alternate.
- *Arms*: rho-neg (reference), ic-iso-rootfind, ic-iso-mitm, ic-kummer-rootfind (p only), ic-plain, rho-neg-aa.
- *Runs*: requested 1,080 (720 measured); completed 1,080; verified 1,080; failures 0.
- *Audits* (`ecbench verify --replay-all`): 8 × `OK`, 720 replays reproduced, 0 problems (`audits/*.verify.log`).

#table(columns: 3, stroke: 0.4pt,
  [*session*], [*id*], [*audit receipt SHA-256*],
  [p20], [`ECBS1hf8c2188cc07d`], [`a4f34d14930f68df9816f7010ef965b167a06465ecdd1e667c4555da4ad89d45`],
  [p22], [`ECBS1h05103c4dcbe7`], [`827897c3e3ada1525822f8db15ed0eba5ef26326a8974407ac1cf181b244256b`],
  [p24], [`ECBS1h197f65094d83`], [`afc37c4e35439df22524b83797e5abc0c11790637a0e6956848781d8e46bf6df`],
  [p26], [`ECBS1h90f9a3ae74a7`], [`c44a6d99320d5b4006d39282e823df5e1644c8067ad48761078792428f7766a5`],
  [p28], [`ECBS1he5bc53b935d6`], [`7820cd0a8fc147c46831771fbab6905fc6f0ef22e4d74c5647e4fd1377eb3628`],
  [q20], [`ECBS1h27d00df033ec`], [`71cd7c0360e54381c0f30c83dc2746b652a8423ea40800fc7637e091143e723e`],
  [q22], [`ECBS1h61cfb60ca261`], [`591e5a990f5e7c38e1121568b47b2299503761e62ce06f43b87c5c1f55dd4e9b`],
  [q24], [`ECBS1h4f6f61e97a73`], [`1966ffc7794ab8d463f0bc870cc49e3088cb270459fd6ba7c737421a3f4e3c21`],
)

Measured values (unit $S$ = GAE / $sqrt(r)$, 16 measured runs per arm) are in `tables/*.table.md` and `tables/summary-arms.tsv`; paired ratios with 95 % cluster-bootstrap intervals in `tables/*.compare.*.txt` and `sessions/*/comparisons/`. Every IC figure is a lower bound (unpriced native work). Wall times in the tables are L0 and not admitted here. Panel 3: 352 `ecbench fb` builds, 0 failed (`tables/panel3.tsv`).

#figure(image("../../figures/s_by_arm.svg", width: 92%), caption: [Mean $S$ per arm and curve, recorded values (sessions above).])
