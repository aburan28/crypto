#set page(paper: "a4", margin: (x: 1.6cm, y: 1.6cm))
#set text(size: 9.5pt)
#import "../../tables.typ": run-section
#align(center)[#text(size: 14pt, weight: "bold")[Run report: N83 precomputed-table rho control, run 1 (38855595549)]]

Workflow `.github/workflows/n83-precomputed-rho-control.yml`, commit `2f7d8c02`, GitHub Actions run #link("https://github.com/aburan28/crypto/actions/runs/38855595549")[38855595549], started 2026-10-16T05:24 UTC. Two jobs (k3300, k600), each on one hosted `ubuntu-24.04` runner (4 vCPU, shared, L0). Arms per job: IC (frozen producer, one process, 16 targets), empty-table fold rho (16 processes), precomputed rho (a) 650 s and (b) = IC max RSS; IC summaries re-analysed in run 2 with the fixed analyzer. Targets: the 16 public points of `gha-38597929540-1`. Requested 16 targets per arm per job; completed 16 per arm per job, no failures or caps. Raw records: `records/run-38855595549/`. Protocol: `PROTOCOL.md`. Observations only; the verdict is in `REPORT.pdf`.

#run-section("run 1 (38855595549)", "records/run-38855595549")

#figure(image("../../figures/online-ms-run1.svg", width: 100%), caption: [run 1 (38855595549): online interval per target; diamond = median with 95% bootstrap interval.])
