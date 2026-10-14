# Greedy factor base through the `ic` tool, end to end

Run 2026-10-14 on the release `ic` binary built from this branch. Every target below
is a synthetic known-answer or public-hash target, so the evidence scope is the tool's
`synthetic_known_answer` / `exact_yield_census_with_synthetic_validation`. Nothing
here is a claim about a real target, and no rho comparison was run.

## What `ic` now does with the greedy family

| use case | command | covered by |
|:--|:--|:--|
| pick a greedy base for a curve, validate it, save the recipe | `ic search --family greedy … --spec-out R.json` | new `Family::Greedy`, `--family greedy` |
| solve a known-answer target with a saved recipe | `ic run --factor-base R.json …` | recipe kinds `greedy_coverage`, `explicit_orbits` |
| precompute the factor-base log database | `ic logs --factor-base R.json --database DB.json …` | same recipe loader |
| recover a target's logarithm by descent from the database | `ic solve --logs DB.json …` | unchanged |
| staged, resumable select → collect → logs → solve | `ic workflow --params wf.json --dir DIR` | `"family": "greedy"` in the params |
| equal-column comparison against random orbits and the shipped sampler | `cargo run --release --example koblitz_greedy_factor_base_yield` | the example |

`--family greedy` is opt-in and is not part of `--family all`, so existing searches
score exactly as before.

**Recipes are frozen.** The greedy scans the whole prime-order subgroup, which costs
seconds on a degree-17 curve. A recipe that stored only the algorithm's inputs would
pay that on every `run`, every worker and every resume. `ic search` and the workflow's
select stage therefore save `explicit_orbits` (the chosen orbit abscissae, one per
orbit), wrapped in `pruned` where pruning applied. Before this, validation runs with the
greedy recipe took 2.2 s each on K₁/2¹⁷ (re-running the scan); after, 0.02 s, the same as
the shipped sampler. A frozen recipe also keeps meaning the same base if the selection
rule is later improved.

## Results

### K₀/2¹³ (`r = 2003`, exhaustive targets)

`ic search --degree 13 --curve-a 0 --family greedy --union-min-seed 3 --union-max-seed 4 --exhaustive-cap 70000 --validate-top 2 --holdout 2 --no-saturate --spec-out k0n13.json`

- selected: `pruned` over `explicit_orbits`, 3 columns, 78 points, coverage 0.896, expected trials 6.7
- validated on 2 holdout fixtures each (known-answer solves), median process time 0.011 s
- `ic run --factor-base k0n13.json --solver pair-table --known-log 777`: recovered 777, verified,
  3 columns, 78 points, 3 relations, 4 trials, 0.010 s
- the legacy `--factor-index 0` base on the same target: 77 columns, 4005 points, 16 trials, 0.34 s.
  That base is a weak baseline, so this gap is not the greedy's effect; the controlled
  comparison is in `REPORT.pdf`.

### K₁/2¹⁷ (`r = 65587`, exhaustive targets)

`ic search --degree 17 --curve-a 1 --family greedy --union-min-seed 2 --union-max-seed 5 --exhaustive-cap 70000 --validate-top 2 --holdout 2 --no-saturate --spec-out k1n17.json`

- best expected trials 18.8 (13 columns, 442 points, coverage 0.850). The same search with
  `--family subgroup` (the shipped sampler's candidates, same budgets and seed) scored 19.8,
  about 5% more. Pruning narrows the gap the controlled experiment shows between unpruned bases.
- `ic run --factor-base k1n17.json --solver pair-table --known-log 4242`: recovered 4242, verified;
  13 columns, 442 points, 14 relations, 16 trials, 0.009 s
- `ic logs … --database db17.json`: complete, 13 columns, 59 relations, 64 trials
- `ic solve --logs db17.json`: `--known-log 40001` and `--random-target --seed 11` both recovered
  and verified. `--known-log 99991` is rejected as it exceeds the subgroup order, which is the correct
  behaviour.

### Staged workflow, K₁/2¹⁷ (`wf17.json`, `wf17-report.json`)

Five targets: two `known_log`, one `random_seed`, two `public_hash_seed`.

- all four stages complete; all 5 targets recovered and verified
- select 6.44 s (the search scoring greedy candidates), collect 4 ms, logs 6 ms, solve 5 ms
- the persisted `factor_base.json` is the frozen recipe; resuming the finished run takes 24 ms

## Limits

- Two curves only (K₀/2¹³ and K₁/2¹⁷), synthetic targets, pair-table oracle, pair sums.
- The greedy still enumerates the subgroup (`r ≤ 2²²`), so `--family greedy` yields no
  candidates on larger subgroups; the search then simply has fewer candidates.
- `--family greedy` scans every orbit up to `r ≤ 2¹⁷` and a sample of 64 per step above that;
  the sampled form loses part of the gain (REPORT.pdf, construction cost section).
- The select stage of a search pays for building every greedy candidate once (seconds at
  degree 17). That cost is paid once per curve and recipe.
