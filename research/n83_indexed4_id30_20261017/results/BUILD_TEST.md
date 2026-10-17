# N83 four-byte ID30 release build and exact-key test

The opt-in four-byte ID30 root-index unit test passed in an optimized
release build on an ARM64 macOS 26.6 host. It checked the largest encoded
state ID, nonzero slot encoding, first-witness retention, and full canonical
key comparison after slot and one-bit fingerprint collisions. The built
executable SHA-256 is
`b3fd21350dc988e3f8980fab870439a1169f2204a990200dba514503c6f69849`.
The frozen K4096 comparison plan SHA-256 is
`14584ea80b9e1a93fc5d3ef0d9ce6de712ae1628b20e78819ecfbf310da1690b`;
its fresh-point selection plan SHA-256 is
`44a068cefe27ae4d9e51220f836bcdccf7d1d86bc480b8219abbc70f0101a275`.

The [build/test receipt](build_test_receipt.json) records the exact Cargo
commands, toolchain, source and binary hashes, and hashes of the versioned
[test stdout](cargo_test.stdout), [test stderr](cargo_test.stderr), and
[build stderr](cargo_build.stderr). The versioned test stdout removes only
one trailing blank line; the receipt also binds the exact raw stdout
SHA-256 `8fd2d0e5ba0d16f9c0faf6e29bac993f2697ead6716442e10daec9ace2d92163`.
The test returned one pass, zero
failures. The source's Lean proof map and kernel/axiom audit are retained
in the parent experiment. The K4096 memory and target-online comparison
has not been run; its sealed 128-point screen and fresh 24-point holdout
remain separate gates.
