# IC improvement round 3: results (2026-10-15)

Two follow-ups to [`../ic_improve_20261015/README.md`](../ic_improve_20261015/README.md), run under [`RULES.md`](RULES.md).

- Each candidate is measured against the incumbent component it replaces.
- Comparisons are paired, with training and held-out splits.
- Every logarithm or root set is verified.

| # | Candidate vs incumbent | Verdict | Key numbers | Landed |
|---|---|---|---|---|
| [01](01-genus3-relations/README.md) | Genus-3 IC relation phase vs round-2 item 06 `sched:lpsgebw2`. **R1:** cheaper per-trial primitives, with identical relations. **R2:** R1 plus a re-tuned κ. | R1 TRADES; R2 ADVANCES (constant) | <ul><li>Whole solve: 0.635× (R1) and 0.617× (R2) at 2^22.</li><li>R1's relation phase: 0.471×.</li><li>Exponent: R2 matches the incumbent (Δ +0.0008).</li><li>Local exponent at 21–22: 1.376 (R1) vs 1.409.</li><li>The instrument is widened to 2^22.</li><li>192/192 logs verified.</li></ul> | example only |
| [02](02-harvest-gate/README.md) | Koblitz Gröbner linear-row harvest: root-rank gate (`gpv < 0.9607` or root rank ≥ 331) vs round-2 item 05's gpv gate | **ADVANCES** (constant), no regression | <ul><li>1.41× at (29,15) m2, 1.89× at (37,19), 1.18× at (31,16).</li><li>m3: 1.04–1.05×.</li><li>End to end: 1.20× and 1.14×, with verified logs.</li><li>0/248 root-set failures.</li></ul> | **Default** in `koblitz_groebner.rs`. `KIC_LINEAR_HARVEST=gpv` restores the previous gate; `=0` turns harvesting off. |

## What the genus-3 diagnosis changed

The relation phase is not what keeps the exponent above 4/3.

- **Trials grow below 4/3.** At fixed κ, the number of trials grows as q^1.258.
- **The log factor is small.** The per-trial log factor adds about 0.05 to the exponent.
- **The κ schedule drove the earlier figure.** The 1.43 seen in round 2 came from the falling κ schedule.
- **The schedule falls because of linear algebra.** κ drops because linear algebra grows faster than q^{4/3}. Roughly, the exponent at the optimal κ is 4/3 + (2/3)(per-trial log growth) + (1/3)(LA excess).
- **Next lever:** the linear-algebra exponent, currently about q^1.52–1.54 locally.

## Caveats

- **Host:** a shared 4-vCPU container at isolation level L0.
- **Measured unit:** operation counts; wall-clock times are descriptive only.
- **Scale:** sizes are small. 2^23 projects to just over the 60-minute job cap.
- **Conductor:** the control plane was unreachable.
