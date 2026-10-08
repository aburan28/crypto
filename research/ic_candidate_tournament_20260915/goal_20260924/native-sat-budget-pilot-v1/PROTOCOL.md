# Disclosed SAT conflict-budget pilot, version 1

This is a fixed **two-query diagnostic**, not a new natural-yield panel or an
attempt to repair the consumed 512-query CryptoMiniSat registration. The
original ordinary F5/CMS target-free capsules and all three historical
confirmation sets remain closed. The accepted CMS executable is invoked
separately on read-only, archived source instances; the original execution
files and audit results are never modified or replayed as a worker.

Hypothesis: raising the native CMS cap from 100,000 to 1,000,000 conflicts,
with a 120-second per-role watchdog, may resolve a selected feasible query and
a selected proved-negative query that were both inconclusive at the original
cap. Selection is explicitly based on the **observed** F5 mathematical
classification, so these two points cannot estimate population yield, rank,
cost per row or any speedup. Run both file and prestarted-stdin modes through
the existing source-validating `ic_sat_stdin_probe`, in the fixed order below.
Preserve timeout, budget exhaustion, model rejection, parity failure and raw
stdout/stderr exactly. Do not retry a failed cell under this protocol.

| Order | Original ordinary query | Public point | F5 geometry | Required CMS outcome for a pass | Export manifest SHA-256 | CNF SHA-256 |
| --- | ---: | --- | --- | --- | --- | --- |
| 1 | 000 | `[103925,114545]` | proved negative | `SOURCE_UNSAT` | `c5efa26c543d0c0efddea60f0103472af3cdafcf48a86f98bae04710e37aee05` | `ddde57d4e0b355d3f01773f128efe7ebda54defcdd81b6a59ab3729546ed4baa` |
| 2 | 015 | `[103877,110181]` | verified feasible | `SAT_MODEL` with source-valid ANF/CNF and full-point lift | `c52941652ba9981a5c59a5470f5a35245122c8944e385370e312dd0e0474c465` | `0e337d0a9b15579c087a94b61fd81940fa5e42b3009fc5e7156038a013e1f3a1` |

The CMS binary must hash to
`6c509f09622f103d8a3ad90afc151e1c4031275c052d7f8af481465b4f27f2af`.
The probe source SHA-256 is
`7dba0b049fb225d074450411d42b0eada200fa43972cd78f3d2be5d4d7050f50`;
the separate geometry preparation is
`67319e5bb0c65e9ae9a10bfbd929487afd5a7e4e8c885b7d123da79b385a9d49`.
Pin the built probe binary before the first role. The probe verifies CMS,
manifest, ANF/CNF/Magma bytes, native source syntax, model satisfaction and
full-point lifting, and retains both original role outputs. It drains child
process groups but its current result explicitly says
`native_child_group_drain_audited=false`; do not promote that to a complete
source-bound execution or target boundary. Run under the shared `busy`
serializer on this physical Mac. Timing remains exploratory.

If either query remains inconclusive, that is evidence against the proposed
budget for this point, not a license to increase it on the same cell. Any
future full-rank SAT candidate needs a separately frozen source/configuration,
fresh natural panel and independent audit; the old F5 rows cannot fill its
matrix. A SAT one-target controller, same-point paired references and host
isolation still follow that gate.
