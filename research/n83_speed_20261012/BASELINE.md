# Baseline and scope

Baseline: main commit c62555b7a1746490dba707589d791669a20bb0d0,
which includes the prior tuned-control study from PR1722. The follow-up branch is
codex/n83-rho-speed-20261012. The unchanged production StrongRhoG implementation
is used by both monomorphizations in examples/n83_rho_speed_eval.rs; the scalar
candidate is isolated in this new example. No earlier measured example, binary,
fixture, raw record or PDF is changed.

The historical retained optimized binary is available at
/workspace/n83_optimization/frozen_binaries/rho_optimization, SHA-256
21827dcc645fd4591b39dd913619f047d267a98d541854bd7c83ee7143d8bbcb.
Its selected configuration is Fold64/DP6/escape-off. This campaign remeasures
that same algorithm and configuration with reference u64 scalars in the new
executable on fresh fixtures; historical times are not a denominator.

Conductor checks were attempted by root and implementers before editing and
returned command-not-found. No .conductor workflow/task files exist in this
checkout. Work ownership is recorded by source directory: new scalar example,
native_runner, verification, audit, reports, plus root protocol/evidence index
and benchmark index. No external messages or deployments are authorized.

Research code, fixture generation, monitoring, selection, statistics and
verification are native Rust/C++. Python is limited to rendering/plotting
already recorded values, PDF inspection, and evidence-archive packaging.

The exact legacy a=1 representation is retained by user context. It is not the
a=0 primary ECC2K confidence-gate fixture. Common calibrated S, an ecbench L2
admission, independent-host replay and challenge scaling remain outside this
measurement. The generic counting scale and automorphism quotient are unchanged.

Initial requirements: a faster exact degree-83 rho implementation, a fair
comparison with the previously optimized rho control, fresh held-out evaluation,
complete costs and correctness, and an updated comprehensive PDF. PROTOCOL.md
freezes the bounded experiment and its acceptance conditions before timing.
