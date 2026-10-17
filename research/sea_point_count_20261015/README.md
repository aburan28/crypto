# SEA and point counting review, 2026-10-15

Requested: find algorithmic improvements in the existing SEA and point-counting
implementations; determine whether CM curve generation exists and implement it if absent.

The implementation is in `aburan28/crypto/isogeny_algos`. Baseline revision:
`cf2745214d83e7bbf737f07d758e486d6af8fa7f`. The repository was located through the
connected GitHub account because the initial workspace had no checkout. No other
repository was modified.

| Requirement | Status | Evidence |
| --- | --- | --- |
| Find algorithmic improvements | Implemented and verified by relevant correctness tests | Adaptive recombination, sorted intersections, invariant scalar reuse; [report](reports/REPORT.md) |
| Check CM generation | Present before this change | `cm::hilbert_class_poly`, `gls_loop::{find_instance, build_curve}` |
| Implement CM if absent | Condition does not apply | Existing constructor retained; independent exhaustive-count regression added |
| Measure an end-to-end speedup | Not performed; no speedup claimed | Correctness validation only; proposed larger changes remain unmeasured |

See the [source-linked report](reports/REPORT.md), [diagram](reports/flow.svg),
[PDF](reports/report.pdf), and [validation log](validation.txt).

Reproduce the relevant checks from the repository root:

```sh
cargo test --manifest-path isogeny_algos/Cargo.toml --release --lib --test sea --test api --test kernel
cargo clippy --manifest-path isogeny_algos/Cargo.toml --release --lib -- -D warnings
rustfmt --edition 2021 --check isogeny_algos/src/find/sea.rs isogeny_algos/src/curve.rs isogeny_algos/src/gls_loop.rs isogeny_algos/tests/sea.rs
python3 research/sea_point_count_20261015/reports/render.py
```

The Python file only renders this recorded report and its diagram; all curve
arithmetic, searches, and correctness checks execute in Rust.

Final validation: 40 tests passed, one pre-existing ignored library test,
Clippy with warnings denied passed, formatting passed, and the PDF's three
pages were rendered and visually inspected. Runtime and performance crossovers
were not measured.

