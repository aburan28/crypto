# Offline freeze validation protocol

Question: can the new target-free controller produce an internally consistent
offline source/dependency/build capsule, with compiled source/build pins and the
accepted external CMS archive, before any scientific invocation?

This is a build and source-custody control on the disclosed educational n17
implementation, not an experiment estimating solver yield, timings or speedup.
The committed `build-validation-config.json` uses three planned queries only
to exercise the bounded configuration schema. **Execute none of them.** Invoke
`ordinary-control-freeze --validation-only` from committed source, then the
frozen worker's `build-identity` command only. The generated validation seal is
not a scientific preregistration. Old target controls and confirmation sets
remain closed.

Use the existing busy lock for bootstrap and offline release builds. Build the
development controller with incremental compilation and debug information off
to bound disposable storage; the freezer clears its child environment and
records its own release arguments and environment. Put the fresh capsule under
`/private/tmp` on the separate system volume, using a new attempt directory.
Use the accepted Homebrew Cargo/Rust compiler explicitly. Freeze inventories
the copied source and vendored dependencies before/after build, verifies the
worker's compiled identity, all build receipts, complete immutable tree and
independently extracted native asset files. No exporter or CMS executable runs.

Success requires the native freezer to complete, all compiled pins to match the
recorded source/build descriptors, and the full source/build/native-asset checks
to pass. Stop at any error, retain the original output and partial capsule, and
fix its cause before a separately named new validation attempt. Do not resume
or overwrite a failed attempt. A failed build is not a solver failure or yield
result. Disk/resource limits and failures stay explicit.

Retain configuration, host context, source commit, lockfile, registration/seal,
identity and build logs/receipts, immutable inventory and a hashed archive of
the immutable tree. No build target/cache is scientific evidence. The archive
is data custody only; never launch its archived worker as a scientific run.
Document the exact commands and any omitted build cache separately. No timer
from this build answers the single-target IC question.

Actual F5/CMS scientific registrations, their fixed 512-query preparation
panels, frozen execution audits, new fresh target/IC1/workload admission,
independent one-target verification and matched incumbent/rho comparison remain
separate gates after source review. CPU speedup additionally requires host
isolation and noise evidence; this control supplies none.
