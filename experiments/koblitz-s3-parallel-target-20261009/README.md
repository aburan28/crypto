# Cooperative target search for the S3 point-decomposition index

Fourteen concurrent, rotating starts in the prepared S3 index reduced the
N53 one-target online wall time by an exploratory **2.67× geometric mean**
against the prelifted group-sum path on the same 24 frozen public points.
Fourteen of the 24 point-level paired ratios were at least 2×. Every run
recovered and verified its scalar; the checked Sage build independently
recovered all 192 run scalars and checked 714 distinct relation witnesses.
The host A/A noise gate failed, so the controlled wall-time ratio remains
unknown. The N41 one-worker control had a 0.925× exploratory ratio.

## Exact candidates and workload

Both curves use `y² + xy = x³ + 1` in polynomial basis, a natural low-`x`
scanned factor base, signed Frobenius folding into 244 effective columns,
and 243 target-independent prepared relation rows. The `fb` counts are actual
distinct subgroup-usable points before folding. Each process receives one
previously unseen public point. The scalar used to construct it is outside
the timed interval. The [frozen plan](FREEZE.json) identifies each workload,
input law, seed, source and binary hash. The [parent plan](../koblitz-s3-crosscurve-20261009/FREEZE.json)
holds the 200 point pool per curve; T001–T024 per curve were selected before
measurement.

| Curve | Field modulus | Subgroup order `r` | Actual `B` / columns | Group-lift candidate | Parallel candidate | Online target workers |
| --- | --- | ---: | ---: | --- | --- | ---: |
| `EC1N41Ckb1hcbf4a37f6e66` | `x^41+x^3+1` | 549756390943 | 20008 / 244 | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0hd7645060ef00` | `IC1N41Ckb1fb20008PDP4rootRCsampleLAgaussTDdirectISO0h8a0e7bb5a56a` | 1 |
| `EC1N53Ckb1h738a572b91f2` | `x^53+x^6+x^2+x+1` | 21044858204113 | 25864 / 244 | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0hcc3dbaef9275` | `IC1N53Ckb1fb25864PDP4rootRCsampleLAgaussTDdirectISO0he884cf865574` | 14 |

The full candidate manifests in `identities/candidates/` retain generators,
base digests, exact solver and code identities, and the per-curve worker policy.
The old and new arms each use at most 14 relation workers during reusable
preparation. The new N53 target query starts 14 workers at spaced index
positions; the first verified relation wins, cancellation stops further
search, and every worker is joined before the one-target clock stops. The
reported S3 root-call count sums every worker. N41 uses the same one-worker
query sequence as group-lift and provides a source/control check.

## Paired local panel

For each point, two runs of each arm follow ABBA or BAAB order. The online
clock starts after factor-base, index, prepared rows and input loading. It
includes target validation, query, all concurrent point-decomposition work,
relation checks, target recovery and scalar replay. Thread startup and drain
are charged. The runner verifies the exclusive phase sum, point, scalar,
rank and factor-base size, and retains failed runs. Local hardware was
Darwin arm64, Mac16,7, 14 physical/logical CPUs, 48 GiB RAM, Rust 1.93.1;
host-wide exclusive partition, fixed frequency and IRQ routing were not
attested. [Hardware record](local/hardware.json).

| Curve | Native / Sage verified runs | Group-lift / parallel median online ms | Exploratory paired geometric mean ratio (target bootstrap 95%) | Points at least 2× (Wilson 95%) | Group-lift / parallel median target S3 calls | Group-lift / parallel A/A p95 abs log ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| N41 | 96 / 96 | 0.484 / 0.507 | 0.925× (0.731–1.135) | 1/24 (0.007–0.202) | 982 / 982 | 1.149 / 0.706 |
| N53 | 96 / 96 | 14.214 / 4.174 | 2.672× (1.731–4.055) | 14/24 (0.388–0.755) | 26569 / 14407 | 0.911 / 2.036 |

The 10% A/A gate is `abs(log(repeat1/repeat2)) <= log(1.10) = 0.0953`
at p95. All four measured p95 values exceed it; the host also lacks the
required isolation receipt. The local ratios are algorithm diagnostics, not
controlled CPU wall-time speedups. Setup medians were 404.1/396.3 ms for
N41 and 2989.8/2912.9 ms for N53 (group-lift/parallel), outside the online
clock. Peak RSS was unavailable in these local rows. The
[analysis](local/analysis.json) retains all 48 per-point paired blocks,
run IDs, failures, timing, root-call ratios and uncertainty.
Across the 24 shared target labels, N53's arm effect divided by N41's is
2.887× (paired-label bootstrap 95% 1.827–4.505). The worker policy itself
differs by curve (14 versus one), so this estimates the chosen per-curve
implementations rather than an isolated field-size effect.

Using the N53 local between-target log-ratio variance, a two-sided 95%
interval with 80% power for a geometric mean above 2× gives a **provisional
112 point-pair** estimate. The 24-point target bootstrap interval still
includes 2×. The high A/A noise makes that sample-size estimate unsuitable
as a fixed Modal run count: start with the frozen 24-point pilot, apply the
host fidelity gate, and recompute the requirement from its paired variance.
N41's current estimated effect is below 2×, so this estimator gives no
finite sample count for that threshold.

## Reproduce and inspect

The parallel source SHA-256 is
`71f527cf7d4f401ad252a6f5ad3fa504c3f8becac9567e254e10d0b48cd75cfa`;
its local release binary SHA-256 is
`ff090babc4477419dac6bac051670fd75e1e3832b34da1886839042a5d7c87c1`.
The plan SHA-256 is
`5733cc4709cbfac6387ccd7514c5ddef19342a4aa5e0c0ab22aaae097f29b5a9`.
The compressed [raw ledger](local/raw_paired24.jsonl.gz) contains all 192
reports and failed-run slots; its uncompressed SHA-256 is
`d11a4c77444b8d07afab825af8bd4bddb24cae458130b44ad6fe26ac1b4bea09`.

From the repository root, with a release build under `CARGO_TARGET_DIR`, run:

```sh
python3 experiments/koblitz-s3-parallel-target-20261009/run_paired.py \
  --first 1 --count 24 \
  --binaries /private/tmp/s3-prepared-target-20261009/release/examples \
  --output /private/tmp/s3-parallel-reproduction
python3 experiments/koblitz-s3-parallel-target-20261009/analyze_pilot.py \
  --ledger /private/tmp/s3-parallel-reproduction/raw.jsonl \
  --output /private/tmp/s3-parallel-reproduction/analysis.json
```

The checked Sage launcher was used for [independent replay](local/sage_replay.json),
with its [runtime receipt](local/sage_runtime_info.json) saved before replay.
The replay verifies each target point against the declared generator,
each target relation against the factor-base points, every distinct prepared
relation equation, the matrix rank, and the independently recovered scalar.
All 24 N41 points retained the same target relation across arms; only two
N53 points did, as expected from independent rotating search starts.

The next measurement is the frozen N41/N53 Modal panel with CPU topology,
affinity, NUMA, PSI, throttling, interrupts, steal time, memory and optional
`perf` IPC telemetry. Guest diagnostics can score noise, while promotion of a
CPU wall-time claim requires the host-level isolation record specified in
the repository rules. The [Modal panel](modal_parallel_panel.py) builds both
sources on the VM and runs the same 24-point ABBA/BAAB design. Its local
`--preview-upload` inventories 1,264 files (52,415,856 bytes): Cargo inputs,
selected `src` and `docs` trees, two S3 examples, three frozen plans, and
48 public target-point files. The preview and generated remote code compile
locally; a parser fixture built from actual T001 reports passed through the
[Modal analysis](analyze_modal_parallel.py). The cloud execution has not yet
been authorized.

```sh
python3 experiments/koblitz-s3-parallel-target-20261009/modal_parallel_panel.py \
  --preview-upload
# After explicit approval to transfer the inventoried files:
python3 experiments/koblitz-s3-parallel-target-20261009/modal_parallel_panel.py \
  --first 1 --count 24 --output /private/tmp/s3-parallel-modal-panel.json.gz
python3 experiments/koblitz-s3-parallel-target-20261009/analyze_modal_parallel.py \
  --inputs /private/tmp/s3-parallel-modal-panel.json.gz \
  --output /private/tmp/s3-parallel-modal-analysis.json
```
