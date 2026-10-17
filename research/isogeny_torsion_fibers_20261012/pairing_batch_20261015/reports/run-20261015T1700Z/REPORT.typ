#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Pairing batch V1 failure and V2 charged retry
*Cell:* $(p,r)=(61,19)$. *Host:* shared Darwin arm64, L0. *Inputs:* `RUN_SPEC_V1.json` and `AMENDMENT_V2.md`.

V1 wrote its setup row and exited on the first unverified rho seed. V2 retained that seed and charged retries. Target zero's first attempt ran 806,452 steps and 386,302,417 current-thread CPU ns without verifying; retry one solved in 11 steps and 32,292 ns. V2 then completed all 18 pairing targets. Its rho total is dominated by that failed attempt.

#image("GRAPH.svg", width: 100%)

V1 raw SHA-256: `175c1ad24d0d43090a3fe5e737880c2f5b5450bf74153a98cb28cb615a084017`. V2 raw SHA-256: `20cc208f41a2dda7d3ee6dae1822d64eebde91d56a06ccb8ccae09884bb9e000`. Both stderr files are retained. The V3 amendment froze an 80-step cap for this subgroup and a separate six-class reference before new runs. These V1/V2 records do not support a method cost ratio. See `REPORT.md` for links and limitations.
