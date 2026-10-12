# Large-degree isogeny GPU component (bounded prototype)

This component constructs a cyclic odd-prime-degree Vélu isogeny when a
**rational kernel generator is available**, then evaluates the map for a batch
of input points on CUDA. It is a native C++/CUDA prototype over prime fields
`5 <= p < 2^31`. A kernel point is found in the CPU fixture by counting the
curve and clearing the cofactor. The GPU accepts an already certified kernel;
it does not discover a kernel, an arbitrary target, or a path to that target.

## Contract

- Source is a nonsingular short Weierstrass curve over a prime field.
- Odd prime `ell != p`; `G` has exact order `ell`, and every nonzero
  multiple of `G` is supplied once. The caller validates these conditions.
- The GPU evaluates `x(P) + sum_Q (x(P+Q)-x(Q))`, and the corresponding
  `y` formula, for each input independently. Infinity and kernel points
  map to infinity. The host computes the Vélu codomain from kernel sums.
- `p < 2^31` bounds every unreduced product in `uint64_t`, and the field
  primitives perform modular arithmetic without overflow in that range.
  The implementation is variable time, intended for public research inputs.
- Cost is `O(batch_size * ell * log p)` field operations in the GPU kernel
  (inversions are exponentiations). No square-root Vélu or compact
  representation is implemented. A prime degree of 101 is a correctness
  fixture, not a performance claim or a P-256 experiment.

## Build and verify

```
make -C gpu/large_isogeny test
make -C gpu/large_isogeny gpu-test ARCH=sm_90
```

The CPU test searches nonsingular curves at `p=1009` for a degree-101
rational kernel, checks its exact order and distinct points, computes the
codomain, checks its point count, and checks image membership and the
homomorphism identities on 30 source points. The CUDA test compares 1024
images with the CPU implementation and checks codomain membership. The
CUDA test must actually run before reporting GPU correctness or throughput.

## Requested large-prime P-256 route

The established Rust `isogeny_walk` already computes small-degree P-256
neighbors through modular polynomials, Elkies kernel recovery, and an
independent kernel verifier (`docs/curves/ic/README.md`). This prototype does
not lift that pipeline to large `ell`. P-256 has prime `#E(F_p)`, so it does
not supply rational order-101 kernel points for this component. Its large
prime-degree horizontal kernels generally require Galois-stable subgroup
construction from modular polynomials or another ideal-action representation.

Next implementation obligations:

1. Specify whether the desired input is `(E, ell)` (enumerate horizontal
   neighbors) or `(E, E')` (find a connecting path). These are distinct
   problems with distinct search costs.
2. Replace the current dense `Phi_ell` and kernel algorithms with algorithms
   appropriate for the target degree, retaining independent certificates.
3. Add 256-bit field arithmetic and an extension-field/kernel-polynomial
   representation to the CUDA path. This point-list kernel is too large and
   often not rational over the base field.
4. Measure cold construction, failed candidates, kernel/map verification,
   transfer, and GPU evaluation separately on frozen degree and batch grids.
   No ECDLP improvement is established by constructing an isogeny.

References: [Sutherland, *Isogeny volcanoes*](https://msp.org/obs/2013/1-1/obs-v1-n1-p25-s.pdf),
[Galbraith, *Climbing and descending tall isogeny volcanos*](https://doi.org/10.1007/s10959-024-01378-4),
and the existing `src/cryptanalysis/isogeny_walk/` implementation.

## Evidence (2026-10-07 UTC)

| Requirement | Status | Evidence / gap |
|---|---|---|
| Native prime-degree map construction | Verified in the rational-kernel fixture | `make test`: `p=1009 degree=101 source=(1,7) order=1010 target=(753,118) generator=(573,570) verified_points=30` |
| CUDA batch map evaluation | Implemented, unverified | `batch.cu` and `test_gpu.cu`; this environment has neither `nvcc` nor an accessible GPU |
| Large-degree P-256 isogeny discovery | Open | No large-degree modular-polynomial/ideal search, 256-bit kernel representation, or P-256 certificate |
| GPU acceleration measurement | Open | No CUDA runtime or matched timing here |

The change adds no new verified curve node or route; the canonical curve
graphs and rendered reports therefore have no new graphable result.
