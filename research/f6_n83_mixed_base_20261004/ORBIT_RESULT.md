# Exact Frobenius closure of the n83 K0 terminal base

The [preregistered orbit screen](ORBIT_PROTOCOL.md) closed the same
cofactor-projected dimension-12 K0 base under all 83 Frobenius powers and
point negation. Every seed orbit returned to its seed. The pinned generator
satisfied `π(G) = [λ]G`, and the exact set consisted of 2,001 complete
signed orbits, each of size 166. The native process exited 0; its
[raw row](orbit_raw.jsonl) and [build log](orbit_run.stderr) are retained.

| Base/set | Actual usable subgroup points | Signed Frobenius columns | Four-summand count / `r` | Five-summand count / `r` |
| --- | ---: | ---: | ---: | ---: |
| Original projected dimension-12 base | 4,054 | 2,027 sign-only columns | 4.662×10⁻¹² | not used for this gate |
| Frobenius-and-sign closure | **332,166** | **2,001** | 0.000209791390 | **13.937281219707** |

The closure has 81.935× as many usable points as the seed base and
0.987× as many relation columns. It passes the registered five-summand
**capacity** gate, after capping the count/order ratio at one. These are
combinatorial upper bounds for a uniform subgroup target. Repeated group
sums can reduce coverage arbitrarily; no ordinary relation yield or
decomposition cost was measured. The dimension-12 base's 2,027 columns
fold only point negation because that set is not Frobenius-closed. The
closure's 2,001 columns fold sign and all Frobenius powers, whose log
coefficients are fixed by the proved subgroup eigenvalue. Cofactor
projection commutes with Frobenius, so orbit-closing the projected base
does not introduce points outside the subgroup.

One **post-screen design calculation** is to keep the existing 8,219,485
entry terminal pair index from the 4,054-point seed base and take the
other three summands from the 332,166-point closure. Its multiset-count
ceiling is
`C(332168,3) × C(4055,2) / r = 0.020765058572` (about 2.08%).
This was derived after the registered all-orbit gate and is only a
capacity bound. For a *fixed choice of three Frobenius charts*, the
three-summand S4 residual system has `3×12+83 = 119` Boolean variables,
which fits a 128-variable monomial representation. The full search still
has chart choices, a terminal-pair support constraint, and many possible
residuals; no existing solver joins them without potentially enormous
enumeration. The 119-variable block is an architectural lead, not a
working n83 F6 solver.

The native command was:

```sh
RAYON_NUM_THREADS=1 CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target timeout 900s cargo run --offline --release --example f6_wide_n83_k0_orbit_capacity > research/f6_n83_mixed_base_20261004/orbit_raw.jsonl 2> research/f6_n83_mixed_base_20261004/orbit_run.stderr
```

Source SHA-256:
`93a0228b6107af8e3eb14dca1d115ee173812b4be9604bb0dae0f23f24743d35`.
Protocol SHA-256:
`bf49b36c9da92f51fd3f134efab898c4e5b9d5e81bdb2d13549313c0817285c9`.
Raw JSONL SHA-256:
`84f2d3421e883cd09d004d17f2412c438c7c89f366cb3d007606adcc7330d477`.
The projected closure set BLAKE3 digest is
`db3e7f25877c75dfd6f53db46447972ce6a44434cb7ea1e53e450d21af0fe79c`.
Its encoding is the same sorted 33-byte affine encoding defined in
[the mixed-base result](RESULT.md). The run was unisolated and contains
no CPU speed claim.

The next admissible step is a source-pinned compositional solver that can
intersect exact three-summand orbit support with the terminal pair index
without scanning all terminal sums per chart. It must distinguish found,
refuted, exhausted, and unsupported outcomes, replay every relation in
the curve group, and first pass exhaustive small-curve controls. The first
n83 performance gate remains one ordinary verified relation under a
frozen budget. A complete one-target IC/rho result is still unknown.
