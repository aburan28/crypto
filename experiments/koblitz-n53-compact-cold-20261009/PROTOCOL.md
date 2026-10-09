# N53 compact-orbit cold full-rank comparison

The first cell will build a fresh 244-column signed-Frobenius orbit base, solve
every base logarithm by guided four-summand S3-root relations, and recover one
held-out public point. The paired strong signed-Frobenius rho run must recover
the scalar of that exact point. This tests whether the compact orbit index
removes the n53 setup bottleneck observed in the signed-expanded producer.

## Frozen mathematical input

Use the binary curve `y^2 + xy = x^3 + 1` over the polynomial basis
`GF(2)[x]/(x^53 + x^6 + x^2 + x + 1)`, prime subgroup order
`21044858204113`, cofactor `428`, and generator
`[198217578752339, 7929897206038174]`. The one primary target is
`[6322155974735900, 5110849281364254]`. Its verification scalar is
`7948768810114`, obtained from SHA-256 of the ASCII string
`n53-compact-cold-20261009|primary-target|v1`, reduced modulo `r-1` and
increased by one. `verify_workload.py` reconstructs the point with independent
polynomial-basis arithmetic; the strong rho preflight also produced the same
point and scalar. Only `target_points.jsonl` is passed to the IC executable.
The scalar is a replay sidecar and an input to rho's point construction, which
is excluded from both online intervals.

The candidate is `construct:53:0:244`, with rank seed `20261009` and no
isogeny. Require exactly 25,864 distinct usable base points and 244 orbit
columns from the emitted base receipt; a mismatch invalidates the proposed
configuration. The compact ID remains unresolved until the executed source,
exact base digest, and complete candidate manifest are frozen.

## Source and execution gate

The current `crypto` default branch has a duplicate rho CLI backend guard
that rejects `strong` before dispatch, and the release library build has
unrelated duplicate Rust definitions. Run measurements only after the
cleanup and one-line CLI fix are merged, the exact source commit and lockfile
are recorded, both examples compile, and their example tests pass. Freeze the
source and executable hashes in a receipt before collecting timings. This
protocol fixes the input and decision rules; it is not a timing receipt.

`run.py` implements two separate steps. After the final source builds, run
`freeze --ic-binary <absolute path> --rho-binary <absolute path>` and save its
JSON output as `freeze.json`; commit and push that file to this PR **before**
running either arm. The freeze records the source commit, hashes of both
executables, the lockfile, the two source files, this runner, and all workload
inputs. Then invoke `run --freeze <absolute freeze.json> --ic-binary <same
executable> --rho-binary <same executable> --run-dir <new run directory>`.
The runner rejects modified or uncommitted source, changed binaries or inputs,
and a freeze file that differs from the committed version. It refuses an
existing run directory, records a preflight attempt, and preserves each arm's
stdout, stderr, exit/resource status, and hashes. It always runs the frozen
rho point after the IC arm, including when IC fails. The 16-GiB gate polls
`/bin/ps` and also records the child's `wait4` peak. Before either arm begins,
the runner checks that `/bin/ps` is usable in the actual execution context;
a denial is a recorded preflight failure. This sandbox denied that check in
ordinary execution, while an approved monitor preflight succeeded.

For IC, invoke `koblitz_orbit_dlp_fast_online` with
`construct:53:0:244 target_points.jsonl 20261009 <output.jsonl>` and set
`KIC_DUMP_BASE` and `KIC_DUMP_RANK` to unique output paths. For rho, invoke
`koblitz_rho_fixture 53 0 signed_frobenius 1 strong 20261009 7948768810114`.
Use default strong parameters (32 lockstep lanes, four distinguished bits)
and one process per arm. Record all effective `KIC_*` variables, CPU model,
thread count, executable and input SHA-256, RSS, stderr, exit status, and
whole-process wall time. Clear algorithm-affecting ambient overrides. Run IC
first and rho second. Cap each arm at 900 seconds and 16 GiB RSS; preserve a
timeout, OOM, or failed verification as its own row with no replacement seed.

The primary measure is verified one-target online wall time. IC starts at
`target_query_begin` after base, index, and log setup and stops at
`recovery_check_end`; sum its five exclusive target phases. Rho starts at its
first target-dependent walk computation and stops after scalar replay. Pair
only when both report the frozen Q and independently verified scalar. Record
IC curve/base/index/rank/matrix/LA costs as a supplementary cold total, with
all rank attempts and failed decompositions charged. Check the cold phase sum
against the in-process interval. A missing phase or failed replay leaves the
end-to-end comparison unknown. CPU time on an unisolated host is diagnostic;
a promoted timing ratio requires the repository's isolated-host receipt.

After a successful run, invoke `replay.py --run-dir <run directory>
--workload workload.json --require-receipt --out <new verification.json>`.
It independently replays every emitted rank relation, modular row and rank
transition, every base orbit, the target group relation, and both final scalar
multiplications. It also checks hashes of all raw files named in the receipt.
Retain the raw relation trace, base dump, JSONL, resource receipts, and replay.
The verifier passed a four-column n13 correctness smoke and rejected deliberate
rank-row and target-scalar changes; that smoke is not an n53 measurement.
After this primary cell is frozen and verified, write a separate protocol for
the secondary multi-target amortization question; do not substitute that
batch result for the one-target comparison.
