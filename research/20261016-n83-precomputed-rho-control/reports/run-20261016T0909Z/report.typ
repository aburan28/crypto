#set page(paper: "a4", margin: (x: 1.6cm, y: 1.6cm))
#set text(size: 9.5pt)
#import "../../tables.typ": run-section
#align(center)[#text(size: 14pt, weight: "bold")[Run report: N83 precomputed-table rho control, run 2 (38878565324)]]

Workflow `.github/workflows/n83-precomputed-rho-control.yml`, commit `44f63a66`, GitHub Actions run #link("https://github.com/aburan28/crypto/actions/runs/38878565324")[38878565324], started 2026-10-16T09:09 UTC. Two jobs (k3300, k600), each on one hosted `ubuntu-24.04` runner (4 vCPU, shared, L0). Arms per job: IC (frozen producer, one process, 16 targets), empty-table fold rho (16 processes), precomputed rho (a) 650 s, (a′) = IC wall, (b) = IC max RSS. Targets: the 16 public points of `gha-38597929540-1`. Requested 16 targets per arm per job; completed 16 per arm per job, no failures or caps. Raw records: `records/run-38878565324/`. Protocol: `PROTOCOL.md`. Observations only; the verdict is in `REPORT.pdf`.

#run-section("run 2 (38878565324)", "records/run-38878565324")

#figure(image("../../figures/online-ms-run2.svg", width: 100%), caption: [run 2 (38878565324): online interval per target; diamond = median with 95% bootstrap interval.])
