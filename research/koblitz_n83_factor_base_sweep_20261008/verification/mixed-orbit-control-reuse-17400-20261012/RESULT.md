# Known-answer solver-access control after the K=600 search: source only

The K=600 natural trial and the earlier K=1,182 trial each exhausted a
5,000-conflict SAT call without a model. The cheapest next diagnostic is one
known-answer 3-small/2-residual control on the already pinned K=1,182/K=2,048
pair. Its [fixture](../mixed-orbit-known-witness-source-20261012/fixture.json)
has an independently checked five-point relation. A valid solver candidate
would establish access to this constructed witness, not natural relation
yield. Another conflict cap would leave even this access question unresolved.

The [new supervisor](reuse_control_supervisor.py) uses the **existing Linux
worker binary**, built and charged in the [K=600 gate](../k600-mixed-orbit-first-2400-20261012/RESULT.md),
instead of charging for a second build. It pins the K=600 audit SHA-256
`387d393b9805f35ee1ca6328e8602133a1b9bc036088b1be1a3ba5a29bb8469c`,
its build receipt SHA-256
`cd302908d34d52412e330fb987e605e288af69d5bf4b50ebbea968e7a2112804`,
the ELF SHA-256
`26dcb2b49388fcb89853eabbda01db71d148620cd465bd5458099e0606dcf606`,
source commit `f3d27f68b6`, the exact known-answer fixture, and both retained
S3 archive hashes. The old 15,000-second source-only control gate remains
unchanged as a historical record; its authorization no longer matches the
current ledger.

The new authorization schema is `n83.additional-local-budget/v5`, purpose
`one_known_witness_control_reuse_k600_worker`, with
`prior_total_charged_seconds="17400"`, the exact K=600 audit and binary hashes,
a fresh output path, an explicit approval reference and exactly **1,200 fresh
local seconds**. The proposed conservative charges are 120 seconds for
read-only `verify-control`, 900 for one `control-run` SAT call, and 180 for
cleanup/audit. A fully launched gate would bring the cumulative local ledger
to 18,600 seconds; a failed earlier stage would be charged only for stages
actually launched plus cleanup. The worker wall is capped at 840 seconds,
with one CPU, 4 GiB, zero swap, read-only inputs and no network. Every
container exit, receipt, input hash and cleanup state is checked; any returned
candidate is independently replayed against the archived point domains and
group equation. A cap, wrong fixture, stale grant or changed evidence cannot
be promoted to a control success.

The proposed control uses the K=1,182 object at
`s3://crypto-autoresearcher/factor-bases/icv1/etc/koblitz_n83_factor_base_sweep_20261008/v2-size-frontier/a0/objects/76e1cbe4ea5a9d26e5b24a2d27f9401271dc471846529190a0416d8e09837a92.jsonl.gz`
and the K=2,048 envelope already recorded in the K=600 gate. The supervisor
requires local downloads with their pinned compressed SHA-256 values before
any container launch; it writes output outside a clean committed checkout.

Three source-only tests passed: the exact new grant is required, a stale grant
cannot reach the container launcher or create output, and the retained audit,
build and binary digests match the committed evidence. The CLI help check
passed. Preparing this gate did not launch Docker, import a retained base or
run SAT; the 17,400-second ledger is unchanged. It cannot select a factor
base or provide a total-runtime comparison. A future run needs a separate
explicit authorization and its own immutable result audit.
