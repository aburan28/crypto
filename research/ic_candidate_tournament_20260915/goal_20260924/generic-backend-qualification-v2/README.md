# Fresh F4/F5/SAT complete-solve qualification

Status: **producer/preservation failure, no campaign artifact** in
[Actions run 36580669479](https://github.com/aburan28/crypto/actions/runs/36580669479).
See the [failure record](RUN-2026092902-FAILURE.md). The seed is closed and
cannot be redispatched; no candidate result or speedup is admitted.
The scientific question, schedule,
resources, success/stop rules and accounting are frozen in [PROTOCOL.md](PROTOCOL.md).
The subsequent [source feasibility audit](STATIC-FEASIBILITY.md) identifies a
hard F4/F5 encoder limit in this registration; it is not a measured campaign
verdict and leaves the live run untouched.
The first registration was [operationally censored](../generic-backend-qualification/RESULT.md);
it established no solver result and must never be redispatched. This panel
uses seed `2026092902`, one process for each of 25 distinct public targets,
250 paired trial slots, and the same five toy curve cells and source-bound
algorithm arms. It does not use the sealed improvement confirmation sets.

The first run's entire possible 25-point fixture schedule is in
[lost-campaign-exposures.json](lost-campaign-exposures.json). An independent full-prepare replay matches that corpus byte-for-byte on all 25 points; see [INDEPENDENT-RECONSTRUCTION.md](INDEPENDENT-RECONSTRUCTION.md). It was reproduced
from the pinned worker and sealed historical exclusions using
[reconstruct_exposures.py](reconstruct_exposures.py); five retained prior
fixtures matched exactly and a second local reconstruction matched the checked
file's SHA-256 byte for byte. The fresh runner passes that corpus to the
tournament's exact curve/point exclusion checker before making new targets.

The workflow is
[`ic-generic-backend-qualification-v2.yml`](../../../../.github/workflows/ic-generic-backend-qualification-v2.yml).
Its PR jobs run the scientific controls and prove that a synthetic interrupted
trial can be packed and uploaded. The original campaign job is now disabled
to prevent reopening this seed. The measured step hit its 300-minute cap;
the original packer deleted the archive after a concurrent-directory `tar`
warning, so the expected campaign bundle was never uploaded. The amended
packer retains and marks an archive from that failure class for future,
separately registered campaigns. It cannot recover this run's lost artifact.

For a future campaign, verify the archive SHA-256, extract it, inspect the frozen
`tournament/evaluator/tournament.py verify` receipt, `natural-yield.json`,
`family-gate.json`, target history, source/build record, and all failure rows.
Run the separately committed [independent_pairs.py](independent_pairs.py)
against the extracted bundle to recalculate the 400 one-target online pairs
from raw receipts and native worker outputs. It checks phase and instruction
closure, target/host/resource pairing, and the published online table. It
withholds aggregate speedups for an arm missing any scheduled paired point.
This post-registration cross-check does not replace the frozen group replay,
natural-query audit, or family gate.
Publish that campaign's bundle link, exact hash and results in a new
evidence PR. A complete, independently verified F4/F5 and SAT arm is the
qualification goal; the development panel is not a held-out or global speed
claim.
