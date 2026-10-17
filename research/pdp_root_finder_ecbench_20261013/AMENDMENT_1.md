# Amendment 1: the oracle's accounting was corrected after the first size ran, and the ladder restarted

Status: written after the first binary's `n = 13` session completed and
its `n = 19` session had started, and **before** any session of the
restarted ladder ran.  The specs, sizes, seeds, arms and decision rules
of `PROTOCOL.md` are unchanged.

## What the first session showed

With the first binary (`cc8c93f8b`), `ic-rootfind` charged 1.26×
`ic-mitm`'s `S` at `n = 13` (41.4 against 32.8, 40 verified runs each,
audit `OK`), against P1's prediction of equality, and its records carried
`lower_bound: yes` where `ic-mitm`'s carried `no`.  Both effects were the
oracle's own accounting, not the measurement:

1. `rootfind-mitm` computed all `2A` points `R ∓ P_k` of a target before
   probing any, so a target that decomposes paid for the whole scan where
   `mitm:m=3` stops at the first hit.  At `n = 13` more than a third of the
   targets decompose.
2. The batched field inversion of each scan was counted as uncharged native
   work (`inversions`), which marks the record a lower bound.  The
   inversion is how the charged additions are computed, not work beside
   them: the cost was complete and was being labelled incomplete.

## The change

Before the restart the oracle computes `R ∓ P_k` in blocks of 32 abscissae
and stops at the first hit (a refutation still costs `2A` additions, a hit
at most one block more than `mitm`'s position of the hit), and counts no
inversion apart from the additions.  `ic_framework::rootfind` and
`docs/ecbench/README.md` §7 say so.  No committed session used the first
binary's `rootfind-*` ids, so the ids keep their names; the first binary's
sessions are kept under `superseded-binary-1/` (the complete `n = 13`
session with its audit, comparisons and tables; the interrupted `n = 19`
session) as the record of why the ladder restarted, and enter no table.

## Disclosure

The first session's other readings were seen before this amendment:
`ic-rootfind / rho-strong = 2.85 [2.64, 3.04]` in `S`, `ic-rootfind /
ic-f4 = 2.56` in `S` (the F4 arm a lower bound) and `0.0030` in admitted
wall, `rho-strong-aa / rho-strong = 1.000` in `S` and `1.016` in wall.
They bear on P2 and P3 and were consistent with them; no decision rule was
changed.
