# N83 K4096 exact one-root local feasibility diagnostic

This bounded diagnostic tests whether the merged IC producer with SHA-256
`3ee6bc11399b382e72789b152043e5a67abc5b783c4e998ad8c534d60baebc03`
can build and solve one exposed N83 point with K4096 folded columns and the
dynamic five-byte, one-root exact index. The public point is the already
exposed cold-screen target with SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
The producer receives only that point; the verifier reads the answer sidecar
afterward and independently replays the factor base, four-point relation and
scalar.

The host is a shared macOS workstation with 48 GiB installed RAM. Require
21 GiB available before launch, stop the child at 21 GiB RSS or below 5 GiB
available host RAM, and cap it at 900 seconds. Use a fresh process and empty
base, index and log system. Save raw stdout/stderr, source and binary hashes,
wall and sampled peak RSS, stop reason and replay. The sample can miss a brief
RSS peak; this run is an algorithmic feasibility diagnostic only and cannot
select a Tier 1 candidate or establish a performance speedup. Preserve a
headroom stop or timeout as such.

The packed four- and five-byte slots encode `state_index+1` in 29 bits.
The source checks this before allocation. K4096 has 696,256,512 symmetric
states, above the 536,870,912-state encoding limit. The exact final K that
fits this state-ID field is 3596; the checked structural proof is in
`../n83_online_range_20261016/OnlineRangeCapacity.lean`.

Run the checked release example built from the current branch with the host's
Python 3.9 environment containing `psutil` 7.2.2. Keep all output on
`/Volumes/SSD990-2/n83-archive-20261016/local-large/`.
