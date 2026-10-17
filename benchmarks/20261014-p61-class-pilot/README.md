# p61 selected-class census pilot

This directory mirrors the visual performance artifacts for the certified
`p61-i0025` selected-class launch study.

Measured on a shared, unisolated macOS 26.6 arm64 host with 8 producer threads:

| Segment | Vertices | Generation | Independent replay |
|---|---:|---:|---:|
| `x = 0..32` | 576 | 13.415 s | 14.778 s |
| `x = 32..64` | 576 | 13.421 s | 14.672 s |
| `x = 64..320` | 4,608 | 21.266 s | 46.939 s c9 / 33.460 s d064 |
| `x = 320..16,704` | 294,912 | 858.472 s | 735.772 s d064; c9 not run |

The measurements are startup-inclusive single observations. A second c9 replay
of the long segment took 215.940 s during concurrent release compiles, so the
timing panel is diagnostic and does not establish a general speedup.
The larger tranche establishes sustained generation and independent replay,
but it has no matched c9 comparison.

The d064 release gate passes the touched binaries, all three continuation-runner
tests, and a complete low-concurrency library run (`4,079` passed, `0` failed,
`114` ignored). Two earlier default-parallel full-library attempts remain
preserved as non-green because their `jv_cover` failure sets varied under host
contention; all exposed failures pass in isolation and in the serial run.

Source of truth:
[`research/p61_isogeny_class_20261014`](../../research/p61_isogeny_class_20261014/).
The full report records the class-order certificate, replay and uniqueness
receipts, storage limits, failures, M7 interpretation, and completion boundary.
