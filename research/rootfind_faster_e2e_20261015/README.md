# Rootfind faster e2e (2026-10-15)

**Decision: engineering** (in progress). Branch
`feature/rootfind-faster-e2e-20261015`. Improves real CPU time of Semaev
root-finding probes and multi-target IC amortization. No security claim;
no ECC2K-130 inference; multi-target secondary to cold single-target.

| | |
|---|---|
| protocol | [`PROTOCOL.md`](PROTOCOL.md) |
| cycles | [`cycles/`](cycles/) |
| host | Apple Silicon; e2e uses `getrusage` CPU time |
| Conductor | T-16 |

## Cycle log

| cycle | change | outcome |
|---|---|---|
| 1 | single-hash probe path; real multi-target e2e | **win**: miss probe ~3×; multi N=8 → 3–5× vs cold |
| 2 | *(pending)* | |
| 3 | *(pending)* | |

## Hard stop vs rho

Cold complete IC still ~1.54× rho-strong (inherited). This branch does not
claim to close that gap; it makes the oracle and amortized multi-target path
honestly cheaper.
