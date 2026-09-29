# Premeasurement preparation failure and amendment

The first local preparation attempt used a sparse checkout that omitted the
source-embedded `research/f4_linear_algebra_20260914/contract.json` file.
`generic_build.py` stopped while hashing source inputs, before compiling or
starting a worker. The sparse checkout was completed without changing its
pinned commit.

The next preparation produced a controlled release build with build-record
SHA-256
`fe6378068da024ffcfdd4d5d23d7d2c1221284ffd6fbe810ca28510177188a93`.
The runner then attempted the first `n17a1/f4` job under the original
`RLIMIT_AS=8 GiB` policy. macOS raised `ValueError: current limit exceeds
maximum limit` in Python's `preexec_fn`; `subprocess.Popen` failed before
`exec` of the worker. Only the 436-byte `job.json` was written. No child
worker, PDP attempt, query, result, or timing observation exists. This is an
operational preparation failure, not an F4/F5 outcome and not a retry of a
measured target.

Before any worker started, the panel was amended to schema v2. The algorithm,
source, build policy, disclosed public points, seed, one-query limit,
60-second wall cap and 8 GiB threshold are unchanged. Instead of an
unsupported macOS address-space limit, the runner samples child RSS every
100 ms and kills a child observed above 8 GiB. RSS sampling is approximate:
a spike between samples can be missed. No pilot result will be used for
performance or memory-comparison claims. The original schema-v1 panel and
runner remain available in commit `691f56580`; the new panel hash is recorded
in the protocol and runner before execution.
