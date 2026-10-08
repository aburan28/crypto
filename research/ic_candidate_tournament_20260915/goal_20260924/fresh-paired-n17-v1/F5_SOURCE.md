# MatrixF5 source gate for the fresh n17 card

The [first target-free publication and independent replay](f5-source-v1/result-v1/RESULT.md)
passed. This is source custody only: card adoption, fresh target execution and
the original post-run audit have not occurred.

The earlier complete F5 control recovered one disclosed point, but its target
capsule was consumed and its full publication included that old target config.
It cannot be reused as a fresh-point registration. This gate builds a **new**
validation-only F5 worker/controller from one clean committed source snapshot,
with the independently audited 512-query F5 preparation bound as a read-only
prerequisite. The validation fixture may be the old disclosed point; it is
never executed by this gate.

`icprog target-control-publish-source` reads the original validation capsule,
verifies its sealed source/build and F5 preparation, and publishes only its
`immutable/` source, vendor, receipts and executable tree. The complete archive
excludes `config.json`, the point-bearing validation fixture, and every target
claim. It retains the validation registration as a small sidecar because that
record contains the immutable inventory and build identity, but no point or
seed. It issues a distinct external **source** seal. Use
`target-control-replay-source` to check every archived member as data without
executing a worker or an archived program. The future `f5` point-card
descriptor must bind that source seal, source manifest and archive hash.

This source gate alone does not authorize a fresh F5 dispatch. After the
four-arm card exists, `target-control-adopt-card` runs from this **same new,
unconsumed validation capsule**. It irreversibly records the adoption start,
preserves its disclosed validation sidecars, and writes a fresh config and
scientific registration for the card while leaving every prepublished
`immutable/` byte in place. `target-control-audit-adoption` must pass before
the scientific registration is published. The F5 worker may then run once
under its original frozen controller and full scientific publication. That
controller replays the source/card adoption before its claim and again in its
original post-run audit. A changed source manifest, binary, card or descriptor
blocks the arm. The old consumed F5 target capsule is never reused.
