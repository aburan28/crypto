# RESEARCH_TPU_IC — a TPU backend for index-calculus relation work

Preregistered protocol for a TPU (JAX/Pallas) implementation of the two
index-calculus stages a systolic array could plausibly help with:
relation collection (pair-table meet-in-the-middle) and the GF(2)
relation algebra. Written to the `boundary, table, ratio` discipline of
[`AGENTS.md`](../../AGENTS.md).

**Status: code landed and CPU-verified correct; no TPU has run.** Every
speed cell in the table below is therefore *unset*, and this document
makes no speedup claim. It is the §8 "stage diagnostic" posture made
explicit before any measurement exists.

**Language note.** This backend is Python (JAX/Pallas) by explicit user
direction (2026-10-01), which overrides `AGENTS.md`'s no-Python rule for
this specific work. The provenance is labelled throughout. The native
Rust+PJRT port remains the candidate permanent artifact if the device
numbers ever justify one.

---

## 1. The boundary, stated before measuring

A TPU result has to be priced against the hardware the repository already
uses, in the repository's unit.

- **Reference (collection).** The existing pair-table relation collector,
  CPU (`src/cryptanalysis/koblitz_index_calculus.rs`) and CUDA
  (`gpu/ecc2k/pairtable.cuh`), in **calibrated operations** and, secondarily,
  in matched-hardware wall time. These are the numbers a TPU must beat to
  earn the word "faster", and they are *already* where this repository's
  relation collection lives.
- **Reference (whole method).** Pollard rho on the same curve, `S =
  total_operations / sqrt(n)` (`AGENTS.md` §2), which is a flat `S ≈ 1.3`.
  A TPU that accelerates one stage changes `S` only through that stage's
  share of the whole, and the residual-walk thread (§5) is the standing
  warning that a stage crossover is not a method crossover.
- **Floor (collection).** The pair table is `|F|(|F|+1)/2` curve additions
  regardless of backend; a backend cannot do fewer. The TPU question is
  only the **cost per addition**, and on a TPU that cost is dominated by
  the bit-matmul-mod-2 reformulation below, not by a CLMUL.

The boundary is **not** "the TPU kernel ran". It is these reference
numbers, measured on the same instances, which this protocol does not yet
have for the device.

## 2. One unit, one table — all device cells pending

Unit: calibrated operations (primary), with `S = ops / sqrt(n)` for any
end-to-end row. Matmul-cell counts and wall time are secondary and only
from a real TPU host.

| variant | stage | correctness | ops / matmul-cells | wall (TPU) | class |
|:--|:--|:--|:--|:--|:--|
| CPU pair-table (Rust) | collection | ✓ (existing) | *reference* | — | baseline |
| CUDA pair-table | collection | ✓ (host-verified) | *reference* | — | baseline |
| **TPU pair-table (this)** | collection | ✓ CPU-oracle | **pending** | **pending** | **unclassified** |
| **TPU GF(2) linalg (this)** | relation algebra | ✓ CPU-oracle | **pending** | **pending** | **unclassified** |
| rho (matched) | whole method | ✓ | `S ≈ 1.3` | — | reference |

"✓ CPU-oracle" = agrees element-by-element with the scalar oracle
(`ic/reference.py`) over the self-check suite; it is a **correctness**
result, not a speed one. The class column stays `unclassified` until a
device measurement lets §3 decide advance / engineering / relabelling /
accounting. Nothing here may be read as an advance.

## 3. Falsification target, declared in advance

The TPU backend is worth pursuing past this correctness milestone **iff**,
on a real TPU host:

1. **Collection:** matched-instance total operations (or matched-hardware
   wall time, 95% paired CI excluding zero) for the pair-table build +
   `m = 3` collection come in **below** the CUDA reference on the same
   curve, factor base, and target count — with every relation verified in
   the group and zero verification failures; **and**
2. that stage's share is large enough that the method's `S` (every phase
   priced, cold) falls — not just the stage's own cost.

It is abandoned if, at the m=83 confidence gate (`AGENTS.md` §8a), the
bit-matmul reformulation's `O(n^2)`-per-multiply overhead is not repaid by
batch throughput, i.e. the TPU collection cost stays above the CUDA
reference at every batch size the device admits. Inadmissible, as always:
moving the unit, dropping a phase, choosing favourable instances, or
quoting a matmul-cell count as a method speed.

## 4. Why these two stages, and the honest fit

A TPU's systolic array multiplies dense matrices. The mapping this backend
rests on:

- **GF(2^m) squaring and reduction are `F_2`-linear** → a fixed 0/1 matrix
  → an integer matmul mod 2. The **Fermat inverse** is a static
  square-and-multiply chain of those matmuls, so the **batch inversion that
  dominates the pair-table build becomes a batched bit-matmul** — the one
  place collection is genuinely array-shaped.
- **GF(2^m) multiplication is `F_2`-bilinear** → one `(a ⊗ b)·T mod 2`
  contraction. This is where the fit is *poor*: a single multiply is
  `O(n^2)` work for what a CLMUL does in a few cycles. Only batch can repay
  it, and whether it does is exactly the open device question (§3).
- **GF(2) relation algebra is an integer matmul mod 2 outright** — the
  array's native shape, and the better-fitting of the two stages.

Explicitly **not** claimed: that a TPU helps the real ECDLP log solve,
which is over `F_r` for a 95–131-bit prime `r`; modular arithmetic mod a
large prime is a poor MXU fit and is out of scope. The GF(2) linalg here
is the parity/Semaev-style flavour, offered because its *shape* fits.

Two stages also keep phases separate (§5): collection and algebra are
priced apart, and neither is allowed to stand in for an end-to-end claim.

## 5. Phases and what dominates

Collection has the pair-table build (`|F|²/2` additions, each a batched
bit-matmul chain for the inverse) and the `m=3` scan (`|targets|·|F|`
subtractions + a sort/join). The sort/join is **not** matmul and is named
as host/vector-unit work, not charged to the array. The algebra stage is
one matmul plus elimination. An exponent fit over ≥4 field sizes, and the
share each phase takes of the method, are **pending the device** — stating
them from CPU interpret-mode timings would be the §6 "extrapolation
presented as a measurement" error.

## 6. What does not count (restated for this thread)

- A matmul-cell count quoted as a relation-collection speed.
- CPU interpret-mode wall time quoted as a TPU speed (it is neither the
  device nor the metric).
- A batch size chosen to flatter the `O(n^2)` multiply.
- Any speed cell above filled in from anything but a frozen device run.

## 7. Scoreboard

No row is added to `docs/index-calculus-scoreboard.html` yet, because there
is no measurement to add — only correctness. The first *device* run that
produces a frozen operation count rides into the scoreboard in the same PR,
with its class set by §3, per `AGENTS.md` §7. Until then the canonical
record is this file plus the green self-check suite.

## 8. Reproduce the correctness gate

```bash
cd tpu
python3 -m pip install -r requirements.txt   # CPU wheels suffice
./run_selftest.sh                            # JAX CPU + Pallas interpret
```

49 checks: field arithmetic vs oracle, curve law vs oracle, pair-table
keys vs an independent oracle build, end-to-end `m=3` relations verified in
the group, GF(2) linalg vs NumPy, Pallas kernel vs `jnp`, and bit-exact
`pack` / `pair_filter_hash` against the Rust constants.
