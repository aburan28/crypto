# SAT target-free source build and exporter control

This is a validation-only source/build result. It contains no public target
card, target query, SAT solve for a target, one-use claim, online interval or
speedup. The original 512-query natural CryptoMiniSat preparation was
rechecked as a prerequisite; those 512 observations are not new measurements
in this result.

The source commit was `7135af53971e0a34b643b15dbb66c65694e8485f`.
The validation registration seal is
`df824d634e829032cb209374d1f395ced43d0604496414cb5d67f464c3e5c865`.
Its source manifest is
`abe85aa0324f179f984c210fa0b2f426b39ab10792a5e6e8b7758581c8358cbd`.
The frozen build copied the committed Rust source and locked offline vendor
tree, built the SAT worker, `icprog`, and one-request prepared exporter, and
pinned the separately audited marked CMS binary. The full
[`capsule.tar.gz`](capsule.tar.gz) archive has SHA-256
`5eb6ff20746bf89cd7c1a29ffb49486641ec63dd845114c124c1954eb63d2d76`.
The [publication descriptor](PUBLICATION.json),
[registration](registration.json), and loose receipt/control sidecars in
`immutable/` allow `icprog sat-target-replay` to verify all **6,974** archive
members as data. A second replay from these repository bytes returned
[`PASS_DATA_ONLY_SAT_TARGET_BUILD_CUSTODY`](independent-replay.json), with
zero archived binaries executed. Appending one byte to a separate copy of the
archive caused replay to reject it as `SAT published archive bytes differ`.

The capsule's prepared exporter hash is
`85e87a5a060a10e7201c9c02337a9a3dcda3fd1cee458d64e65cdf37d69b5b11`.
Using the accepted exporter hash
`4da7f5781da1aba2f76d6b9ce8f6e4ea43845f2465cf8edf0612615931d0629a`,
the three disclosed points `[40991,73355]`, `[73003,104622]`, and
`[59775,2910]` each passed exact ANF/CNF-XOR/Magma byte parity, deterministic
manifest and source-instance parity, pre-input READY, malformed/oversize
request rejection, and never-fed child cancellation. The
[data-only audit](exporter-parity-audit.json) returned
`PASS_DISCLOSED_EXACT_BUILT_EXPORTER_PARITY`. It checked all three original
results, fixed source hashes, raw stdout/stderr and receipt hashes, exact
request digests, and five matched child PID starts/completions per point.
The [raw controls](exporter-controls.tar.gz), SHA-256
`d711c0beab70347b9583125e7d74302b54cc6f94b166e5f1ef3766207f5a96f9`,
are retained separately. They establish transport parity on disclosed
fixtures, not natural relation yield on a fresh point. The control probe's
binary hash was checked, but no separate process-execution attestation for
the probe is claimed; the independent audit replays the retained data.

The source freeze deliberately rejects a scientific registration at this
stage. A scientific freeze must bind parity for its own exact exporter and a
new one-use SAT controller before it can join the four-arm point-card
campaign. The host was an unisolated physical Apple ARM64 machine, so these
checks do not establish a CPU timing ratio.
