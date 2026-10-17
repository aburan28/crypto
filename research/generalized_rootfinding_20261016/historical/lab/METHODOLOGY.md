# A continuing root-finding research program

The objective is a new exact root-finding method over prime and extension fields
that achieves the frozen 2× end-to-end performance gate. An unsuccessful method
is an experiment result, not a reason to end the program. The repository skill
is [.agents/skills/rootfinding-autolab/SKILL.md](../.agents/skills/rootfinding-autolab/SKILL.md).

The active field targets are GF(2^83) and the standard NIST P-192, P-224,
P-256, P-384 and P-521 prime fields/curves. Other prime and extension fields
remain in the protocol scope. Exact NIST field and group proofs are provided
by `lab/research/nist_fields.py`; the prime route uses interval/implicit-map
domains rather than assuming a proper prime-subfield prefix space. Field size,
coordinate admission and a small-domain root check do not establish hardness
or a comparative speed gain.

An agent supplies the scientific work: form a falsifiable hypothesis, derive
coverage and complexity, look for counterexamples, implement a variant, and
review evidence. The local lab coordinator supplies reproducibility: validate
declared source/input hashes, enforce resource limits and dependencies, execute
commands, and preserve attempts. Repeated execution of an unchanged timeout is
not an improvement strategy. The next attempt needs a new mechanism, a repaired
assumption, or a specific diagnostic purpose.

The loop is **hypothesis → derivation → independent oracle → variant → frozen
plan → bounded experiment → review → next queue**. Keep mathematical correctness,
literature novelty and qualified performance as independent assessments. Append
failed and inconclusive runs as well as successes. Never overwrite frozen source,
input selection, or raw evidence to make a result look better.

Use a 15-minute/three-experiment session by default, with 4-GiB, one-core worker
limits and a 90-second initial component experiment. Session limits are budget
controls; they do not relax the [release gate](../performance/GATE.md). Finite
daemon cycles can execute declared plans, but an agent must design new patches
and analyze results. Empty ready queues pause awaiting that work. Missing F6 or
full-width baselines open engineering tasks while independent derivation,
counterexample and codec experiments continue.

Initialize with `python -m lab.runner init --manifest lab/research_program.json --ledger lab/ledger`.
Then use `status`, `validate`, `run-ready` and `resume`
with `--ledger lab/ledger`; `run-ready` requires explicit experiment and wall
budgets. `daemon` has a finite `--cycles` count and `--interval` between cycles.

Prioritize exact quintic support verification, conditional factor search that
reduces the actual surviving dimension, the 512-bit matched baseline port,
prime-curve coverage, and circuit/representation optimization. Details and
promotion rules are in the [research queue](../.agents/skills/rootfinding-autolab/references/RESEARCH_QUEUE.md).

Every promoted algorithm must preserve all original polynomial roots, including
repetitions and nonrational-lift branches, or explicitly declare its narrower
contract. Hard-target success, novelty and ≥2× superiority remain unestablished
until their separate evidence passes. This discipline turns negative results
into actionable constraints without turning them into an impossibility claim.


## Field-specific structural thresholds

These thresholds describe exact source-polynomial representations. They are
not measured time limits, hardness certificates or qualified speed claims.
For an integer interval D=[0,2^d−1], let M=2^d−1. A source polynomial has
signed elementary coefficients bounded by binomial(m,i) M^i. Lift each
coefficient with explicit multiples of the actual coordinate-field prime;
reducing a canonical Fp residue modulo unrelated small primes is unsound.

| NIST field | Maximum d for unique cubic lifts | Maximum d for unique quintic lifts | Proposed critical cubic d |
| --- | ---: | ---: | ---: |
| P-192 | 64 | 38 | 64 |
| P-224 | 74 | 44 | 74 |
| P-256 | 85 | 51 | 85 |
| P-384 | 128 | 76 | 128 |
| P-521 | 173 | 104 | 173 |

The critical cubic family has 3d source bits, close to the coordinate-field
bit length. Its integer-cubic support oracle isolates roots with exact Sturm
arithmetic, retaining multiplicities and avoiding source-interval scans.
NIST odd group-order proofs exclude rational two-torsion; ordinary arity-three
targets outside D admit a separate complete no-ghost proof.

For five sources, heuristic ordinary-target density suggests investigating
d near 40, 46–47, 52–53, 78–79 and 105–106 respectively. These domains can
wrap their fifth coefficient, so shared quotient variables are mandatory.
Fix exact dimensions and targets before timing; retain every holdout case.
The fastest completed matching baseline still supplies the one-second
hardness floor and the full release gate still decides promotion.
