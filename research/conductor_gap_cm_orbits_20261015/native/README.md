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

The adapter's five exclusive online phases include parsing, child teardown,
independent Rust replay and final target verification. Raw backend phase
diagnostics alone exclude that replay. Supplementary cold costs include map,
orbit and pair-table setup. Operation counts still contain unpriced work;
these L0 runs do not support an eligible IC-versus-rho speedup claim.

`folded_recurrence.cpp` implements the generic recurrence
P[k+2] = P[k+1] - [5]P[k], after P[0] and P[1]=[lambda]P[0].
It produces the same scalar orbit on either endpoint without a CM point map.
`folded_d43.cpp` supplies the separate degree-11 CM control.

`run_cold_orbits.sh` mechanically extracts the unchanged arithmetic from
`folded_recurrence.cpp` and compiles `orbit_process.cpp` as a thin process
wrapper. Twenty fresh processes each contain one shuffled pass of all 600
model/arm jobs. The aggregate records preserve every observation; individual
duplicate batch record files are omitted. Process CPU and orbit CPU are
reported separately, with no host-isolation or cache-flush claim.

`evidence-manifest.sha256` at the study root records this checkpoint's tracked
study files and maintained implementation/UI/workflow sources. Verify from
the repository root with `sha256sum -c` and the manifest path. Regenerate after
changing evidence; active constructor outputs are excluded until sealed.
