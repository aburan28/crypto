# N83 K1107/K1112 root-layout cold screen

This screen compares the root11 and root22 dynamic exact indexes at the two
lowest K values from the complete 1102–1120 exposed-point screen. It keeps
the exact GF(2^83) curve `EC1N83Ckb1a1haacb69ee9f37`, IC source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`,
and previously exposed public-point file SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
The usable base size is `166K` and folded column count is `K` at each value.
The layout changes the exact root index, so each layout has its own candidate
identity even when its factor-base points are identical.

Run eight fresh IC processes in symmetric order: K1107/root11 A,
K1107/root22 A, K1112/root11 A, K1112/root22 A, K1112/root22 B,
K1112/root11 B, K1107/root22 B, K1107/root11 B. Run one same-point rho
process afterward. Use rank seed 29, four rank threads, one target thread,
a 480-second cap per IC process and a 120-second rho cap. Require at least
four exposed CPUs, the same 15.5 GiB RSS envelope and 2 GiB free disk.

Retain source/binary hashes, exact target and verifier sidecar, topology,
affinity, cgroup and CPU telemetry, raw producer streams, full base dumps,
relations, phase costs and every failure. Independently replay all base
points, four-point relations and scalars. Require the execution order,
resource envelope, topology and telemetry gates. Each of the four A/B
candidate pairs must have identical base/index identity and cold spread at
most 5% of its mean; both layouts at fixed K must have the same factor-base
digest. A failed, timed-out or noisy arm records `BOUNDED_UNKNOWN` and keeps
its raw receipt.

Rank the four complete pipelines descriptively by the mean full cold
process interval on this exposed point. Keep the promoted candidate null.
The sealed sixteen-point finalist fixture in
`research/n83_e2e_final_holdout_20261016` supplies the independent next
comparison for the leading K/layout candidates. That final protocol will
also keep the primary online and supplementary cold objectives separate.
