# Preregistered relation-row certificate for rotated normal-basis factors

Status: protocol only; no compressed relation-row outcomes have been run.
This is a bounded follow-up to merged [PR #762](https://github.com/aburan28/crypto/pull/762)
(commit `2de218f583d2dddeafdd0180dd02329dc57d53e4`), whose frozen
raw archive SHA-256 is
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`
and receipt SHA-256 is
`8b35cf4ea223b112d987abeed04a0ba8bf04e2f509ec0943eead87bdcdce9eac`.
Its source/input freeze is `FROZEN.json` SHA-256
`77c313cc3db97430e77a797492c32037eecb6189225fc478789d8db6f659edc3`.
The archive is immutable input; this PR does not edit it or the separate
PDP corpus agent's frozen files. The question is whether every saved
full-point witness induces a correct, compact prime-subgroup relation row
with the known Frobenius scalar. No new PDP search, relation rank,
linear solve, target logarithm, or n=131 solver speed is attempted.

## Algebra and exact n=13 roster

Let E:y²+xy=x³+1 over F_(2^n), H=[4]E(F_(2^n)) of prime order q,
τ(x,y)=(x²,y²), and Fi=τ^i(F0). For a witness
`P_i∈Fi`, `sum_i P_i=Q+T` with one of the four rational
`T∈{O,(0,1),(1,0),(1,1)}`, define
`R_i=[4]τ^(n−i)(P_i)∈H` (τ^n=1 on rational points). If
`τ(R)=[λ]R` on H, then

`sum_i [λ^i]R_i = [4]Q`.

The producer and an independent arithmetic verifier must check the
raw `Q+T` witness, τ inverse transport, the per-term
`[4]P_i=[λ^i]R_i` identity, and the final row equation by actual group
law. At n=13, q=2003, λ=89 and the polynomial is `0x201b`; derive the
order `4q=8012` from the Weil recurrence as in #762 and verify
`τ(H)=[89]H` on its frozen generator. No signed or torsion term may be
silently discarded.

Use **all** saved positive coset witnesses from the four *rotated* #762
arms (beta=3/7 × m=5/6), in ascending (beta,m,k,T-index) order. For
every one of the 2,003 H targets in each arm, inspect all four Q+T
entries and preserve misses in a manifest. The archive contains exactly
1,799, 6,307, 1,799 and 8,012 positive full-point sums, respectively:
17,917 distinct supported full-point decisions in total. Fail if the
new row roster omits, duplicates, or invents a positive from #762.
The repeated-F0 controls in #762 are not rotated Fi witnesses and are
not mixed into this λ^i row test.

For each projected nonzero R_i, choose the lexicographically smaller
coordinate tuple of `R_i` and `−R_i` as its canonical log column,
record the sign `s_i∈{+1,−1}`, and add
`s_i λ^i mod q` to that column. Drop only R_i=O, recording its source
index and proving it came from rational torsion; retain a zero row when
all terms cancel. Aggregate repeated canonical columns modulo q, retain
both the unaggregated signed terms and the final sorted nonzero
coefficient list, and compute its group-law point independently of the
raw signed sum. A row passes only if it equals `[4]Q`. The row still
records k, Q, T-index, raw Q+T, all source witness indices and source
points. A signed point cannot be merged with its opposite without
applying its negative coefficient. Count distinct canonical columns,
zero/torsion terms, repeated-column terms, cancelled columns and row
weights (nonzero coefficient count), with full distributions and stage
operation/CPU/wall/RSS costs. No row is called independent or useful
without a separate matrix-rank experiment.

In addition to all #762 witnesses, for every beta/m make two explicit
exception controls outside the saved positive roster: (1) every slot
contains `(0,1)`, giving an all-zero projected row; (2) slots 0 and 1
contain `P∈F0` and `τ(−P)∈F1` for the lexicographically first
non-torsion F0 point, remaining slots `(0,1)`, giving one repeated
canonical column with coefficients `1` and `−λ`. Construct Q and T
from the raw sum by the exact cofactor-four decomposition, then replay
both controls independently. These controls do not enter the 17,917
archive-roster count.

## Optional n=131 planted structural controls

Before any row outcome, freeze a separate SHA-256 input file of at most
four planted tuples for each of `(m,d)=(5,25)` and `(6,21)` under beta=3
in #762's polynomial model `0x800000000000000000000000000002007`.
The schedule domain is
`ECC2K130-ROTATED-ROW-20260925-v1/m/d/tuple/slot/counter`; read digest
big-endian modulo 2^d, take the first nonzero x-mask with
`Tr(x+1/x²)=0`, and lift via the odd-degree half-trace. Use the sign bit
from SHA-256 of the same string suffixed `/sign`. Tuple 0 reuses one
accepted mask across every slot (a repeated compressed column); tuple 1
sets slot 0 to x=0, `(0,1)` and derives all other slots. Tuples 2 and 3
use independent accepted masks. Transport each accepted signed F0 point
to Fi by τ^i. Preserve rejected candidate counters, chosen masks, sign
bits, full source points, synthetic target Q and T in the frozen input
file. Compute `S=sum_i P_i`, `Q=[4^{-1} mod q][4]S∈H`,
`T=S−Q∈E[4]`; check T against the four rational torsion points.
These are planted public synthetic group identities, not solutions to
ECC2K-130 target PDPs. If building this frozen input or independent
replay cannot satisfy a 30-second/512-MiB per-cell budget, preserve a
failed feasibility receipt and keep the accepted result n=13-only. No
post-result mask or sign selection is allowed.

## Freeze, replay, caps and decision

Commit the producer, separate field/group verifier, n=131 input
construction (if admitted), exact source/input SHA-256 manifest, commands,
archive format and caps **before** the first row output. Open this as a
draft protocol PR before that outcome. The producer may reuse #762's
Euclid-inverse group code; independent replay must use #762's separate
bit-serial/Fermat arithmetic (or another separately implemented group
law), validate the #762 archive hash, rebuild all factor lists and
`Q+T` witnesses, then reconstruct every compressed row from source
indices rather than trusting the saved row coefficients. Replay the
explicit exceptional rows and every admitted n=131 planted row.

Bound each n=13 arm by 120 seconds and 512 MiB primary-process peak RSS,
each n=131 cell by 30 seconds and 512 MiB, and independent replay by
600 seconds overall. Enforce wall deadlines in the producer, reject an
over-cap RSS cell retrospectively, and preserve partial outputs and a
quantitative failure receipt. Record field multiplications, squarings,
inversions, point additions, scalar calls, wall, CPU and process peak RSS;
charge source/archive setup, τ back-transport, [4] projection,
canonicalization, coefficient/group-row evaluation, Q+T checks, and
independent replay separately. Host-specific wall/CPU/RSS are diagnostics;
the group identities and exact counts are primary.

Pass only if the frozen roster is complete, every source and compressed
row identity passes in both implementations, torsion and sign controls
pass, and all source/input/archive hashes agree. A failure is a useful
negative result and must remain archived. Passing this gate certifies
row *semantics* and empirical row compression on the frozen witnesses,
not rank, a solver, relation yield on n=131, a Certicom logarithm, or
an end-to-end speed advantage. The next matrix/rank or solver experiment
needs its own preregistered targets and matched costs. Update the
canonical scoreboard and decision ledger with this bounded result.
