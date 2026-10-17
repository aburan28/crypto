# Frozen reusable-pairing-input batch gate (2026-10-15)

Status: preregistered before generating public batches or measuring reuse.
The original requested counts are `N=1,4,16,64` **distinct nonidentity**
targets per group. The [cardinality record](../../../problems/pairing-batch-distinct-target-cap.md)
proves that `N=64` is unavailable on the existing order-19 and order-61
groups. Keep that requirement pending; do not repeat a target to fill a row.

## Fixed cells and batches

Use the two public pairing inputs in
[`INPUTS.jsonl`](../pairing_transfer_20261015/INPUTS.jsonl), with the exact
`(p,r)=(61,19),(271,61)` source curve, G, and first Q. Use the existing
degree-18 and degree-10 evaluation-field constructions and the same H-search
seeds and caps. For each cell, form a public batch with the original Q first.
Generate a `StdRng` permutation of the remaining distinct nonzero scalars
using seed `0x4241544348000000 XOR p`, map them to `[k]G`, and commit **only**
the public Q list as solver input; keep the scalar fixture in a separate
post-solve checker input. The runner must reject any answer field in its
input. Freeze both files and hashes before timed runs.

Run cumulative prefixes `N=1,4,16` in both cells, and `N=60` on `r=61`
as an explicitly additional feasible endpoint. A separately screened cell
with `r>=67`, independently checked nondegenerate pairing and practical
evaluation field must supply the outstanding `N=64` row under its own
versioned protocol. No result on the old groups substitutes for that row.

## Paired work and limits

For each three-round cell, build and certify the field cold, search H once
including all failed attempts, compute `alpha=t_r(G,H)` once, then compute
`beta_i=t_r(Q_i,H)` for every distinct target, solve each finite-field DLP
without the answer, and independently check `[k_i]G=Q_i` in a separate native
process. Report cumulative `C_setup + sum_{i<=N} C_target_i`, memory peak,
and each phase's field/curve primitive counts, CPU and wall. Price N rho
solves on the **same** Q list with seeds
`0x4241544348001000 + 256*cell + 64*round + target_index`.
Keep cold single-target and batch totals separate. Complete cost includes
field construction, irreducibility, H search, alpha, every beta and finite
field solve, failed attempts, and final scalar verification. Stop at the
original per-row `4096` candidate x, `512` acceptable H,
`10^7` extension multiplications, `512 MiB` peak RSS, or `15 min` thread CPU;
retain any partial prefix as censored.

The fixed generic two-class floor is `S=sqrt(pi/4)=0.886226925` operations
per `sqrt(r)` for one target. The reference is same-target rho, summed for
each batch count. Convert all native work before quoting a counted `S`;
otherwise report whole-process current-thread CPU and wall as shared-host
L0 descriptive values. A batch crossover is an amortized, multi-target
observation and does not settle the one-target ECDLP boundary. Preserve the
degree-18 p=61 solver as a separate implementation control until a minimal
degree-9 output-field solver and independent embedding check exist.
