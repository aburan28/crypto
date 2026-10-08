# tpu/rho — batched Pollard rho for prime-field ECDLP (JAX)

Parallel collision search (van Oorschot–Wiener, r-adding walk, distinguished
points) over short-Weierstrass curves `y² = x³ + ax + b` mod p, with the
field arithmetic written so the heavy part of every modular multiply is an
int8 matmul the TPU MXU can run. The prime-field sibling of `tpu/ic/`, held
to the same standard: *no device has run it, and its arithmetic is verified
against a host oracle on the CPU first.*

> **Language.** Python (JAX) by explicit user direction for the `tpu/`
> backend (2026-10-01), overriding `AGENTS.md`'s no-Python rule for this
> work. Labelled in `rho/__init__.py`.

> **No speed claim.** Everything here is a *correctness* result. There are
> no TPU timings. The CPU steps/s the solver prints are a progress meter,
> not a measurement, and must not be quoted as one (`AGENTS.md` §5, §8).

## Why rho and not index calculus

Index calculus does not apply to generic prime-field curve groups; the
ECDLP workhorse is parallel rho. Its parallelism (millions of independent
walks) fits a batch dimension trivially; the per-step cost is one affine
point addition, i.e. a handful of multi-limb modular multiplies plus one
field inversion — and that is where the MXU question lives.

## The one idea

A modmul `a·b mod p` is three big-integer products:

| product | operands | lowering |
|:--|:--|:--|
| `a·b` | walk × walk | batched Toeplitz matvec `[N,2L−1,L] @ [N,L,1]` |
| `q₁·μ` (Barrett) | walk × constant | dense `[N,L+1] @ [L+1,2L+2]` int8 matmul, shared weight |
| `q₃·p` (Barrett, truncated) | walk × constant | dense `[N,L+1] @ [L+1,L+1]` int8 matmul, shared weight |

Field elements are little-endian vectors of **7-bit digits stored as int8**
(0..127 is a *signed* int8, which is what the MXU consumes: int8 × int8 →
int32 accumulate). `L = ⌈bits(p)/7⌉`; secp256k1 is L = 37. The two Barrett
products have a fixed right-hand operand (Toeplitz matrices of
μ = ⌊B^(2L)/p⌋ and of p), so across a batch of N walks they are exactly the
large weight-shared matmuls the systolic array wants. The walk×walk product
has no shared operand; XLA lowers the batched Toeplitz form to a batched
matmul with an inner N=1, which uses the MXU poorly. That is the known weak
spot — see "Open items".

Carry propagation is a `lax.scan` over the limb axis: L sequential vector
ops over the batch, exact for any int32 input including negatives
(arithmetic shift = floor division, so `a−b` normalises to `(a−b) mod B^L`
with the borrow as a `−1` top carry). A fixed number of parallel carry
rounds would be cheaper but is not exact (a run of 128s can ripple the full
width), so the scan is kept. Max intermediate magnitude is
(L+1)·127² ≈ 6.1·10⁵ at L = 37, well inside int32.

Inversion for the affine add is a **product tree** over the whole batch:
`log₂N` batched levels up, one serial Fermat inversion (`p−2` bits, ≈2·|p|
modmuls on a single element) at the root, `log₂N` levels down. 3N modmuls
in total, and the root cost is per step regardless of N.

## Walk

- Table: `r = 32` points `R_j = a_j P + b_j Q`; index `j = x[0] & 31`.
- Distinguished point: `dpBits` low bits of `x >> 7` are zero — the DP
  window is read from limbs 1..3, independent of the index bits
  (`dpBits ≤ 21`).
- A walk that lands on a DP **freezes** for the rest of the chunk; the host
  harvests `(x, c, d)`, records it, and reseeds the walk at a fresh random
  `cP + dQ`. Nothing is lost; the cost is idle lanes for ≤ K steps per DP.
  Keep `K ≪ 2^dpBits`.
- `x == R_j.x` (doubling / inverse case) is flagged `bad`, its denominator
  replaced by 1 so it cannot poison the inversion tree, and the host
  reseeds it without recording.
- Collision `c₁P + d₁Q = ±(c₂P + d₂Q)` → both candidate `k` are tried
  against `Q = kP`.

## Layout

```
tpu/rho/
  field.py   7-bit limbs in int8; Toeplitz product; Barrett via constant matmuls;
             scan carry; Fermat inverse; tree batch inverse
  curve.py   host bigint EC ops; Miller–Rabin; BSGS order in the Hasse interval;
             prime-order curve generator; secp256k1 parameters
  walk.py    jitted K-step chunk: affine r-adding step, batch inversion,
             DP freeze, bad-denominator flag
  solve.py   table build, seeding, DP store, collision → k, verification
  demo.py    end-to-end solve on a generated curve
tpu/tests/test_rho_*.py   CPU self-checks (run by tpu/run_selftest.sh)
```

## What is verified (CPU backend)

1. **Field** — `mulMod`, `addMod`, `subMod` match Python integers at 40,
   64 and 256 bits (secp256k1 p), including edge values; `batchInverse`
   matches `pow(x, -1, p)`.
2. **Walk** — on **secp256k1** (L = 37), the device step equals the host
   reference walk point-for-point and coefficient-for-coefficient, and the
   solver invariant `W = cP + dQ` holds after every step; the DP predicate
   matches the host predicate.
3. **Solver** — prime-order curve generation (BSGS order, primality, Hasse);
   `solveFromCollision` on both the `+` and `−` branches; an end-to-end
   solve on a generated 32-bit curve recovers `k`.
4. Out-of-suite runs (not committed as evidence, correctness only): 40-bit
   curve solved in 1.05·10⁶ steps vs 1.21·10⁶ expected; 44-bit in
   7.9·10⁶ vs 4.8·10⁶ expected. Both recovered `k`, zero bad restarts.

```bash
cd tpu
./run_selftest.sh                 # all tpu/ tests, rho included
python3 -m rho.demo 40 2048 8     # bits, walks, dpBits
```

## Open items

1. **Pallas kernel for the walk×walk product.** The honest MXU formulation
   is an `L×L` VPU outer product reduced along anti-diagonals by a constant
   0/1 matrix `[L², 2L−1]` — shared-weight, but with int16/f32 inputs
   (products up to 2¹⁴). Needs benchmarking against the batched-Toeplitz
   lowering on a device.
2. **Special-form primes.** For secp256k1 (`p = 2²⁵⁶ − 2³² − 977`) replace
   Barrett with the fold `x_hi·(2³² + 977) + x_lo`: one tiny constant matmul
   instead of two.
3. **Root inversion.** Serial and per step. Sliding-window exponentiation
   cuts it ~40 %; on a device, N ≈ 2¹⁶–2¹⁸ walks is what makes the batched
   work dominate.
4. Move the DP harvest on-device (`cumsum` compaction + scatter) so the
   host sees only a small DP buffer per chunk.
5. **ecbench registration.** Before any cross-method claim this must be
   charged in the counted unit the native rho uses (`ecbench-extend`); until
   then it is a stage diagnostic with every speed cell pending.
