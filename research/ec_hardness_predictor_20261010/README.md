# ec_hardness_predictor_20261010 — geometric ML hardness predictor over isogeny classes

Status: **plan frozen 2026-10-10; nothing run.** Result class: level-1
learned-correlation study with level-2 exact certificates; **not an ECDLP
speedup**. `S`, end-to-end cost and speedup are unset.

| file | content |
|:--|:--|
| [`PLAN.md`](PLAN.md) | the end-to-end plan: question and hypotheses, inventory of what exists in `crypto` and `ml-cryptanalysis`, the 14 gaps, label and feature definitions, data-generation tiers and budgets, labelling pipeline, splits and controls, models (tabular + typed-edge GNN), registered decision rules, phased work plan |
| [`DELEGATION.md`](DELEGATION.md) | the self-contained brief handed to the executing agent |

Companions: `research/isogeny_class_difficulty_20261008/` (protocols I-1..I-5,
whose nulls this plan inherits), `research/large_prime_isogeny_degree_20261008/`,
`research/isogeny_walker_engineering_20261008/`, `research/ecbench_yield_sweep_20261004/`,
and the learner at `/Volumes/SSD990-2/ml-cryptanalysis` (`docs/PLAN.md`, milestone M7).
