# Native backend provenance

`group.hpp` and `folded.cpp` start from the existing local CM pilot sources
(`/workspace/prime-gap-100/group.hpp`, `ic_cm/campaign.cpp`), preserved with
content hashes in the historical summary. New code removes target-dependent
seed exclusion, adds the same-crater scalar arm, construction-only repetitions,
work counters and phase diagnostics. `ecbench::cm_folded` checks the executable
hash before launch and uses independent Rust affine arithmetic to replay every
factor-base scalar, retained relation, matrix rank and recovered target.

Build:

```
g++ -std=c++20 -O3 -pthread native/folded.cpp -o target/cm-folded
cargo build --release --bin ecbench --example conductor_gap_study
```

`native/sqrt_velu.cpp` is a separate codomain-only implementation of kernel
power sums. FLINT 3.1.3 product/remainder trees supply fast polynomial arithmetic;
GMP supplies field arithmetic. Link system FLINT/GMP development packages using
`-lflint -lgmpxx -lgmp`. No point-map or dual-map certificate is implied.

`native/search.gp` runs in PARI/GP 2.15.2. `isprime` supplies rigorous primality
proofs. Kernel order and generator order witnesses are explicit; floor-model
certificates are a separate check after construction.

The backend's supplementary phase diagnostics exclude independent Rust replay
and do not yet meet the complete five-phase online receipt contract. Its raw
operation counts contain unpriced work. A successful recovery is therefore
not an eligible IC-versus-rho speedup claim.
