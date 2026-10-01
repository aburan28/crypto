# Qualified measurement plan

New timings follow current `AGENTS.md` section 10. `isolated_run.py` freezes the
protocol, sources and isolation controller, builds and tests under its `busy`
lock, preserves the actual binaries, and then reserves/pins one physical core
for each worker. It records the full isolation receipt, CPU identity/features,
memory, operating system, compiler, source hashes and executable hashes.

An eight-round A/A run precedes A/B for every fixture. Both A/A labels call the
same retained dispatcher. A/B retains all 52 frozen reference solvers, five
original half-word arms, two feature-gated half-word arms and two matched full-word
three-input-XOR controls. There are 61 A/B arms. Model, logical work and assignment
order comparisons remain mandatory.

Failed, refused, contended or missing-isolation runs are sealed. They do not enter
eligible timing aggregates and are not silently retried or pooled. A new campaign
gets a new directory and retains its relation to earlier attempts. The A/A symmetric
97.5-percentile ratio is an additional noise threshold for each comparison group.
No observation inside that spread establishes an improvement.

The local macOS environment can check the arithmetic and feature-gated kernels,
but the repository's Linux affinity/pressure controller cannot qualify timings
there. New local wall-time measurements are therefore disabled. Linux ARM64 CI is
the intended timing host. The repository is public; GitHub documents
`ubuntu-24.04-arm` as a standard ARM64 runner for public repositories, with standard
public-repository use free:
<https://docs.github.com/en/actions/reference/runners/github-hosted-runners>.
This is a runner-availability fact, not evidence that a benchmark has run or that
the selected host supports EOR3. The capability record decides the latter.

The VM's physical neighbours and CPU frequency remain outside the isolation tool's
control. Report that hardware class and the measured A/A spread, and preserve any
refusals or contention. Do not transfer Linux ARM64 timings to Apple silicon, x86,
GPUs or a cryptanalytic pipeline.

Discovery and the 240-fixture full comparison remain separate. The full comparison
requires unchanged timed Rust sources bound to a sealed qualified discovery run,
including the retained 216 fixtures and 24 unused holdouts. A positive full result
still needs unchanged-source confirmation on further unused holdouts. Full-IC,
calibrated-operation and rho costs remain unmeasured throughout this standalone
Boolean experiment.

Schema 3 executes A/A and A/B for one fixture in a single reserved worker and retains the combined stdout plus exact phase slices. n12 uses sixteen repetitions; other sizes use eight. The resource record covers the paired fixture. This adds real calibration/comparison work without padding or changing the 10% resource threshold; see `ISOLATION_ATTEMPTS.md`. The current discovery therefore has 14,640 comparison and 480 A/A observations; the full grid has 146,400 comparison and 4,800 A/A observations. Earlier counts describe their frozen predecessor protocols.
