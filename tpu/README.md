# tpu — JAX/Pallas index-calculus backend

A TPU backend for the two index-calculus stages a systolic array could
plausibly help: **relation collection** (the pair-table meet-in-the-middle
of `gpu/ecc2k/pairtable.cuh`) and the **GF(2) relation algebra**. It is the
TPU sibling of the CUDA kernels in `gpu/`, and it is held to the same
standard those are: *no device has run it, and its arithmetic is verified
against a scalar oracle on the CPU first.*

> **Language.** This directory is Python (JAX/Pallas) by explicit user
> direction (2026-10-01), overriding `AGENTS.md`'s no-Python rule for this
> work. The provenance is labelled in every module. The native Rust+PJRT
> port is the candidate permanent artifact if device numbers justify it.

> **No speed claim.** Everything here is a *correctness* result. There are
> no TPU timings; the stage-diagnostic table in
> [`protocol/RESEARCH_TPU_IC.md`](protocol/RESEARCH_TPU_IC.md) has every
> speed cell left **pending**, on purpose (`AGENTS.md` §5, §8).

## The one idea

A TPU multiplies dense matrices; GF(2) is not obviously its business. The
bridge this backend is built on:

| field / algebra operation | is… | → on the array |
|:--|:--|:--|
| GF(2^m) add | `F_2`-linear | XOR (elementwise) |
| GF(2^m) square, reduce | `F_2`-linear | fixed 0/1 matrix → **int matmul mod 2** |
| GF(2^m) inverse (Fermat) | square-and-multiply chain | a sequence of those matmuls |
| GF(2^m) multiply | `F_2`-**bi**linear | one `(a ⊗ b)·T mod 2` contraction |
| GF(2) relation algebra | — | **int matmul mod 2**, directly |

So the batch inversion that dominates the pair-table build, and the whole
relation-algebra stage, become batched bit-matmuls — the MXU's native
shape. The honest caveat: a *single* GF(2^m) multiply is `O(n²)` here
against a CPU's one-instruction CLMUL, so only batch throughput could
repay it, and whether it does is an open device question.

## Layout

```
tpu/
  ic/
    reference.py   scalar oracle (no JAX): field, curve, pack, hash — the contract
    instances.py   deterministic toy Koblitz curves (enumerable, solvable)
    field.py       batched GF(2^m): add/sqr/mul/inv as bit-matmuls mod 2
    curve.py       batched binary-Koblitz add/double/neg/frobenius/pack
    pairtable.py   STAGE 1: pair-table build + m=3 relation collection
    linalg.py      STAGE 2: GF(2) matmul / elimination / solve / nullspace
    pallas_kernels.py   the bit-matmul-mod-2 Pallas kernel (TPU + interpret)
  tests/           49 CPU self-checks (oracle agreement, end-to-end, contract)
  protocol/RESEARCH_TPU_IC.md   the preregistered, honest protocol
  requirements.txt  run_selftest.sh
```

## Run the self-checks (no TPU needed)

```bash
cd tpu
python3 -m pip install -r requirements.txt   # CPU wheels are enough
./run_selftest.sh
```

JAX runs on CPU and Pallas runs in `interpret` mode, so the identical
kernels are exercised without a device. On a TPU host, pass
`interpret=False` to `bitmatmul_mod2` and install the matching `jax[tpu]`.

## What is verified

1. **Field** — `mul`, `sqr`, `inv` match the scalar oracle for several
   degrees; `a·a⁻¹ = 1`.
2. **Curve** — `add` over all point pairs, `neg`, `frobenius`, `pack`
   match the oracle on enumerated toy groups.
3. **Pair table** — JAX-built keys equal an independent oracle build;
   every stored key is really its pair's sum.
4. **End-to-end collection** — `m=3` relations `R = P_i+P_j+P_k` are
   found and re-verified in the group (nothing trusts the device).
5. **GF(2) algebra** — matmul / rank / solve / nullspace vs NumPy; the
   Pallas kernel vs plain `jnp`.
6. **Contract** — `pack` and `pair_filter_hash` are bit-exact against the
   Rust constants, so a table built here is readable by the Rust/CUDA
   backends.

## What this is not

- Not a TPU *performance* result — see the protocol; all speed cells are
  pending a device.
- Not a help for the real `F_r` log solve (poor MXU fit; out of scope).
- Not a structural-fidelity instance: the toy degrees here are exploratory
  (`AGENTS.md` §8b); the m=83 gate is where transfer toward ECC2K-130 would
  be judged.
