#set page(paper: "a4", margin: 18mm)
#set text(font: "DejaVu Sans", size: 10pt)
#set par(justify: false)
= Native baseline at p = 2^61 − 1
Run record: 2026-10-16. Status: partial requested experiment; completed capped baseline.

The repository's native two-summand relation collector ran with three inverse pairs over the exact requested prime field. All twelve streams exhausted 100,000 target trials with zero relations and no recovered scalar. This is an observation within a finite cap, not a proof of nonexistence or an operation-cost ranking.

#image("execution.svg", width: 100%)

== Inputs and execution
Source library: #raw("de0c9391480155ebdd5a8495d2c79e4508d4795d").
Driver source revision: #raw("c510e49"), pinned in results/source-commit.txt with complete SHA-256 records. The driver calls the unchanged prime_fast backend. The protocol was committed before execution.

Curve: #raw("p=2305843009213693951; a=p-3; b=111; #E=r=2305843010818082053; h=1"). Generator: #raw("(1,1796364306280342617)"). All points are subgroup checked. Exact ICV1 and target encodings are in the raw manifest; registration of the new model is pending.

One synthetic known-answer target is used on each curve. The collector receives only G and Q. Seeds are 2026101600 through 2026101611. Each run stops at 20 relations or 100,000 target trials; the outer process cap is 300 seconds. Both original processes exited zero.

== Verification and control
The separate 20-bit control used 16 inverse pairs, one seed and the same target cap. It collected 33 relations in 62,676 trials; all 33 were replayed by ordinary group operations, and scalar 596767 passed the independent group check. The large-field streams produced no relations to replay. A second large-field execution was byte-identical to the first; this is same-host reproducibility, not independent validation.

Host: Linux x86_64, AMD EPYC 9V74, five visible logical CPUs, 16 GiB container limit. The host capsule is results/host.json. No isolated wall-time figure is reported. Trials and relation rows are stage counters; total operation cost, S, rho reference and speedup remain unset.

== Remaining scope
The screenshot's 1,000/2,000-candidate driver, degree-7 map and original counting definition were not located in the inspected default branch. This new prime-order a=-3 instance is not a transported screenshot instance. The exact candidate experiment remains open in problems/screenshot-factor-base-driver-unlocated.md and its linked backlog card. No large-field solve, candidate ranking or ECC2K-130 transfer is established.

== Sources
Paths are relative to research/prime61_factor_base_20261016: protocol.json; results/control.jsonl; results/prime61.jsonl; results/replay-prime61.jsonl; results/source-and-binary.sha256; results/host.json; README.md. Source code: examples/prime61_factor_base.rs and src/cryptanalysis/prime_fast.rs. The diagram uses recorded counts; unlike-curve workloads are not shown as a comparative graph. Canonical performance dashboards remain unchanged because no complete-cost or accepted-candidate finding is added.
