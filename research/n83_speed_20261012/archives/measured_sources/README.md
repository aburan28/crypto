# Pinned measured-source copies

These files are byte-for-byte copies from source/protocol commit9da9ff2c,
matching the Oct12 measured executable. ARCHIVE_SHA256 binds the copies;
the original compiled-source digests and binary pins remain in provenance.
The Cargo manifest is named Cargo.toml.snapshot so this partial source archive
is not mistaken for an independently buildable crate.

The complete dependency/source tree is the Git revision9da9ff2c, descended from
main c62555b7. Reproduce by checking out that revision in a separate worktree,
using the recorded toolchain, commands, point-only inputs and new output paths.
These selected copies are a browsing/provenance aid, not a full crate closure.

On Oct14, delivery integrated newer main275beb93f at2ea7fd320. Upstream
StrongRho adds diagnostic trail vectors, the one-word Gf2 changes visibility,
and Koblitz IC changes other implementation paths. The measured binaries and
timings were not replaced by binaries built from those later sources. Do not
run the Oct12 path-based frozen_sources.sha256 check against an arbitrary newer
checkout and interpret its expected mismatch as alteration of the frozen run.
Verify it in the pinned worktree; these archived copies give matching content.
