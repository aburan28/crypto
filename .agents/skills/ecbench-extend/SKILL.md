---
name: ecbench-extend
description: Add an ECDLP method, curve construction, factor base or oracle to the ecbench harness so it is charged in the same counted unit as rho, BSGS, kangaroo and index calculus, verified externally, replayable, and stored in the database. Use when a new algorithm or variant must be measured against the existing ones.
---

# Extend ecbench

Read `src/cryptanalysis/ecbench/mod.rs` (module map) and
[`docs/ecbench/README.md`](../../../docs/ecbench/README.md) §7 and §11. Native
Rust only (AGENTS.md).

## A new method

1. **Implement it on `CountedGroup`** (`src/cryptanalysis/ic_boundary.rs`).
   - Every point addition or doubling goes through `g.add` / `g.double` /
     `g.mul` on a `GroupOps` ledger. Keep set-up and search on separate ledgers.
     They become phases.
   - Count everything else (table lookups, canonicalisations, solver steps) in a
     `BTreeMap<String, u64>` under a name ending `_uncharged`. Never drop it.
   - Take the target point and a seed. **Never take the planted scalar.** Return
     the recovered value without verifying it. The runner verifies it in its own
     process.
   - Be deterministic in (workload, seed). No `thread_rng`, no wall-clock budget
     in the search, and no iteration over a `RandomState` map that decides
     anything. The audit's replays must reproduce bit for bit.
   - Put generic algorithms in `src/cryptanalysis/ecbench/generic.rs`, with unit
     tests recovering known logarithms on a prime curve and a Koblitz curve.
2. **Register it** in `methods::registry()`: id `family.variant`, a one-line
   summary stating its expected cost, the code path in `entry`, `applies`, and
   every parameter. A parameter that changes the search gets no default
   (`default: None`).
3. **Dispatch it** in `methods::solve` (or `solve_generic`), returning a
   `SolveReport` with phases, `automorphisms_used`, counters and
   `unpriced_of(&counters)`.
4. **Test end to end:** add an arm to a scratch spec, run it with `--cpus none`,
   and run `verify --replay`. Every run must verify and every replay must be
   `identical`.
5. **Document the charge** in the README §7 table: what is charged and what is
   counted but not charged. If the method trades memory for operations, say so
   beside the table.
6. **Never change what a registered method counts** once a committed session
   uses it. CI replays every committed run, and the change would fail them. That
   is the point: a changed algorithm is a new id (`kangaroo.vow2`), and the old id
   keeps its evidence.
7. **An IC change changes the candidate identity.** A `vs_rho` claim's IC1
   `candidate_id` binds the SHA-256 of `ic_framework/{mod,plugins,solvers,linalg}.rs`,
   `ic_boundary.rs` and `ecbench/methods.rs` at compile time (`claim.rs`,
   `Implementation::this_binary`). Editing any of them gives a new candidate,
   as it should; the pinned identities in `claim.rs`'s tests use placeholder
   hashes and do not move.

## A curve construction

Add a `workload::CurveSpec` variant that calls an existing repository
constructor, and give `call()` its registry `generator_call` form. Check that
`ecbench plan` shows the slug `registered=true` and an EC1 alias. Otherwise
register the curve (`scripts/build_curve_registry.py`, AGENTS.md §11) before
citing it. Word-size limits: `GF(p)` with `p < 2^62`, and `GF(2^m)` with
`m ≤ 62`.

## A factor base or decomposition oracle

Write it as an `ic_framework` plug-in (`FactorBaseBuilder` or
`DecompositionOracle`), price its native work through `price_phase`'s counters,
and add its name to the `match` in `solve_ic_prime` or `solve_ic_binary` and in
`dump_factor_base`. Check that the same plug-in spec gives the same `FB1h…` from
`ecbench fb` and from an `ic.pipeline` run. The database then joins them.

## Before the PR

```bash
python3 tools/isolated_bench.py busy -- cargo test --release --lib ecbench
```

```bash
python3 tools/isolated_bench.py busy -- cargo test --release --test ecbench
```

Run `rustfmt --edition 2021` on changed files and clippy with `-D warnings` (CI
checks both). Run `python3 scripts/check_curve_names.py --diff origin/main`.
