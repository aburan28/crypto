# ecbench smoke: every method (2026-10-12)

Stock harness smokes on the PR tip, run for operation counts on macOS (L0;
no affinity). Not a performance claim.

| Session | Spec | Arms |
|---|---|---|
| `sessions/prime` | `docs/ecbench/specs/smoke-prime.json` | every rho/BSGS/kangaroo/IC method on a 20-bit prime curve |
| `sessions/koblitz` | `docs/ecbench/specs/smoke-koblitz.json` | signed-Frobenius rho, negation rho, BSGS, mitm-frobenius IC on n=17 |
| `sessions/binary-lp` | `docs/ecbench/specs/smoke-binary-large-prime.json` | `ic.large_prime` (n−1, double LP) vs signed-Frobenius rho on n=11 |

```bash
cargo build --release --bin ecbench
./target/release/ecbench run --wait --allow-busy --spec docs/ecbench/specs/smoke-prime.json --out research/ic_techniques_20261012/ecbench_all_20261012/sessions/prime
./target/release/ecbench run --wait --allow-busy --spec docs/ecbench/specs/smoke-koblitz.json --out research/ic_techniques_20261012/ecbench_all_20261012/sessions/koblitz
./target/release/ecbench run --wait --allow-busy --spec docs/ecbench/specs/smoke-binary-large-prime.json --out research/ic_techniques_20261012/ecbench_all_20261012/sessions/binary-lp
./target/release/ecbench verify --dir … --replay 8
./target/release/ecbench table --dir …
```
