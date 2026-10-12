# N41 feasibility smoke, unisolated macOS

This pair preceded the larger-panel protocol freeze. It establishes that the
existing S3 baseline and paired-inversion binaries both accept the N41 Koblitz
instance and recover the same known fixture scalar. Checked Sage independently
replayed the public point and recovered scalar; the 252 relation witnesses were
checked for agreement between the two binaries but were not independently
replayed. It is **not** a pilot, confirmatory sample, or a controlled CPU
speedup result.

The fixture is `123212651130 * G` for `KoblitzCurve::new(0,41)`, where the
constructor's generator is `[2056947637384,1635505394702]` and the prime
subgroup order is `549756390943`. `target_fixture.rs` records the exact
construction; only `public_target.json` was passed to each solver. The source
snapshots are the immutable `baseline.rs` and `candidate.rs` in the adjacent
`koblitz-s3-pair-query-20261007-v4` experiment.

Commands, with the binary paths used on this host:

```sh
/Volumes/SSD990/crypto/target/release/s3-baseline-v4 41 0 244 20260928 \
  /private/tmp/s3_n41_target.json /private/tmp/s3_n41_baseline.jsonl 14
/Volumes/SSD990/crypto/target/release/s3-candidate-v4 41 0 244 20260928 \
  /private/tmp/s3_n41_target.json /private/tmp/s3_n41_candidate.jsonl 14
```

Each process exited 0. The actual base had 20,008 subgroup points and 244
folded columns; both runs reported `group_verified=true`, scalar
`123212651130`, 252 relation attempts, rank 243, and the same target relation
indices and base digest. Their **exploratory** target-online intervals were
181.602292 ms and 82.272583 ms, respectively (baseline/candidate 2.2073).
The process was neither pinned nor isolated and its noise counters were not
captured. The full raw JSON lines and empty stderr files are retained here;
the stdout was the same JSON line as each `.jsonl` file.

SHA-256 of the binaries used:

| Arm | SHA-256 |
| --- | --- |
| Baseline | `145457f4f66a7f54d6b09fac69ad16b2647b1d00ef63f214c9b3de8d0ff92ad4` |
| Candidate | `366e5c08145339aad0b141c689009ce606e9e707f1b0bf96580f17fe4a6f24cf` |

The existing v4 freeze receipt binds those binaries to the exact macOS build
sources. Linux measurements require new binaries and manifests; the macOS
binary hashes are not portable to the Pod.

After the protocol was written, `s3_icms_adapter` was exercised once on this
same smoke point with the frozen baseline binary. The command used the same
`41 0 244 20260928 ... 14` arguments, with
`adapter_metrics.json` as the metrics path. Its raw solver line is
`adapter_raw.jsonl`; the empty `adapter.stderr` is retained. The adapter
reported an exact five-phase online interval and verified scalar 123212651130.
This checks the adapter's translation on one input and is not an additional
isolated benchmark pair.

`replay_n41_point_sage.py` was run through the checked repository Sage launcher
after saving `sage_runtime_info.json`. Its immutable-input hashes and scope are
in `independent_point_replay.json`. The replay confirms the GF(2^41) model,
subgroup membership, fixture target, both recovered scalars, and agreement of
the rank and relation witness fields between arms. It does not validate the
relation witnesses independently or change the exploratory timing status.
