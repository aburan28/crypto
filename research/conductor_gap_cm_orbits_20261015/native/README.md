# Native backend provenance

`group.hpp` and `folded.cpp` start from the existing local CM pilot sources
(`/workspace/prime-gap-100/group.hpp`, `ic_cm/campaign.cpp`), preserved with
content hashes in the historical summary. New code removes target-dependent
seed exclusion, adds the same-crater scalar arm, construction-only repetitions,
work counters and phase diagnostics. `ecbench::cm_folded` checks the executable
hash before launch and uses independent Rust affine arithmetic to replay every
factor-base scalar, retained relation, matrix rank and recovered target.

Build from the repository root:

```
study=research/conductor_gap_cm_orbits_20261015
mkdir -p target
g++ -std=c++20 -O3 -pthread "$study/native/folded.cpp" -o target/cm-folded
cargo build --release --bin ecbench --example conductor_gap_study
```

`native/sqrt_velu.cpp` is a separate codomain-only implementation of kernel
power sums. FLINT 3.1.3 product/remainder trees supply fast polynomial arithmetic;
GMP supplies field arithmetic. Link system FLINT/GMP development packages using
`-lflint -lgmpxx -lgmp`. No point-map or dual-map certificate is implied.
The recorded environment uses FLINT 3.1.3; its internal product/remainder-tree
API is not promised compatible with other releases. The v4 source is separately
versioned, with 4,096 polynomial-coefficient checks and two 100-reference passes.
See [the arithmetic derivation](CONSTRUCTOR.md).

```
g++ -std=c++20 -O3 "$study/native/sqrt_velu_batched_v4.cpp" -lflint -lgmpxx -lgmp -o target/cm-sqrt-velu-v4
target/cm-sqrt-velu-v4 --self-test
target/cm-sqrt-velu-v4 64 < "$study/prepared/kernels-reference.txt"
```

The historical launchers contain explicit paths to this workspace's extracted
dependencies and timed binaries. A run on another host needs a separately
versioned launcher and fresh hashes; do not rewrite old receipts. The guarded
96-bit launcher refuses local compilers or solver tables at launch. Keep them
stopped until the bounded constructor exits; this is resource coordination,
not a physical-host isolation receipt.

Registry updates and website rendering have separate responsibilities. Use the
native `upper-ids` / `index` commands for independently certified model records.
Then rebuild the repository's native cover/trait/transfer sidecars and existing
canonical website views from the repository root:

```
cargo run --release --bin curve_cover_check
cargo run --release --bin curve_traits -- build --budget 262144
cargo run --release --bin jacobian_transfers -- build --budget 200000
python scripts/build_ic_leaderboard.py
node tools/render_ic_dashboard.mjs
python scripts/site/panel_index.py
python scripts/build_lab_browser.py
```

These existing website tools render frozen repository evidence and catalog
joins; new experiment aggregation and statistics remain in native Rust.
The obsolete partial native browser/leaderboard refresh commands were removed
because they could leave existing rows stale when a model gained a generator.
Preserve the ranked board when checking regenerated views. Do not run the
native compilation/sidecar stage concurrently with a large constructor.

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
