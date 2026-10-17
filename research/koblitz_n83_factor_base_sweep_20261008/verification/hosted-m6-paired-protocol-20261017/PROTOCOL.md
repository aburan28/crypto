# Hosted N83 K600 4+2 paired natural-search gate

This protocol freezes one public trial against the same sequential K600 base
and two different K600 pool streams. The Gray pool has 503 shared signed
Frobenius orbits and 97 residual orbits; the hash pool has zero shared and
600 residual orbits. The three compressed archive byte hashes, S3 identities,
public target, trial coefficients, source capacity estimates, worker arguments,
and limits are fixed in [plan.json](plan.json). `plan.py --check` reconstructs
that file from hash-pinned prior receipts and the whole-subgroup support census.

The shape is four base summands plus two residual summands, with the first
pair ordered. Each arm gets the identical public target and trial zero, at
most 20,000 solver conflicts, 64 models, 2,300 worker seconds, and 2,400
external seconds. Archive regeneration and the shared worker build occur on
the same hosted Linux ARM runner using the digest-pinned Rust image. Both
arms run in separate networkless, read-only containers with one CPU, 4 GiB
memory and zero swap. The job has a 7,500-second ceiling. Count that entire
ceiling as a conservative hosted charge until raw job timing is audited.

The worker's `PASS_one_bounded_partial` means a candidate passed its group
equation on this one trial. `UNKNOWN_search_cap`, worker/outer wall caps, and
an absent second arm are censored. `PASS_one_trial_refuted` applies only to
the frozen trial. Producer, import, image, archive, or cleanup failures are
infrastructure/evidence failures. Neither a positive nor a negative one-trial
result estimates relation yield. Total index-calculus runtime and the best
factor base remain unset pending repeated independent queries, useful graph
rank, linear algebra, target log, verification, and matched rho accounting.

The workflow retains the regenerated archives, build logs, worker binary,
container inspection, raw solver outputs, and outer receipts. A separate
post-run replay must independently inspect all hashes, domain overlap,
candidate group equations, solver counters, caps, and S3 round trip before
any empirical statement is promoted.
