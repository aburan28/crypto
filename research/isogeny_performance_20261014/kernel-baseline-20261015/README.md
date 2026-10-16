# Frozen kernel baseline: 2026-10-15

Source `5521605f86ae444a9f4a433ca9b2f80d541e8470`, with native executable in
`bin/`, compiler receipt, exact standalone source archive, two complete raw
runs, stdout/stderr and the native reader's comparison. Both native run shells
exited zero. The source archive builds the performance binary without Cargo
dependencies or external runtime curve data. Cargo remains a maintainer build
tool; users run the captured executable directly.

```sh
RAYON_NUM_THREADS=1 ./bin/isogeny-perf kernels --out FRESH.jsonl \
  --revision 5521605f86ae444a9f4a433ca9b2f80d541e8470 \
  --bits 192,224,256,511 --ells 7,31,101
./bin/isogeny-perf compare baseline.jsonl repeat.jsonl FRESH_COMPARISON.json
```

Each run has 84 cases plus its header and terminal PASS. All 84 paired output
fingerprints agree and all four operation counters repeat exactly. The report
helper verifies that this comparison names the raw run SHA256 digests before
rendering any count. This is a baseline, with zero qualified wall samples and
no algorithm optimization claimed. The generic Field model leaves inversion
internals opaque. Setup, references and verification are excluded; owned result
allocation is included. Different supplied inputs and output materialization
boundaries stay explicit in `../KERNEL_BASELINE_PROTOCOL.md` and the reports.

For each fixture, both Kohel and point-enumerating Vélu references agree on
degree, codomain and full point transports at `2*ell+1` distinct non-kernel
abscissae. Counted and native constructors reproduce that fingerprint. This
check is distinct from the BMSS series' full coefficient comparison. These
are synthetic supersingular field-width fixtures, not named NIST searches.

The source is stable throughout both runs. Actual compiler: Homebrew rustc
1.93.1; host Apple M4 Pro, 14 logical CPUs, 51,539,607,552 memory bytes (48 GiB).
The shared host also had unrelated work and the root validation compile.
No timing samples were requested or promoted. The full standalone native
release suite and strict all-target Clippy pass in the retained logs.

The complete raw directory is in `../KERNEL_BASELINE_RUNS.tar.gz`. Restore with:

```sh
tar -xzf research/isogeny_performance_20261014/KERNEL_BASELINE_RUNS.tar.gz \
  -C research/isogeny_performance_20261014
```

`FILES.json` applies to this extracted directory. Uncompressed local files are
also retained unchanged. Reports and graph data are in
`../kernel-report-20261015/`; matching browsing copies are in
`../../../benchmarks/20261015-isogeny-kernel-baselines/`.
