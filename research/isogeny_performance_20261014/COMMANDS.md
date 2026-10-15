# Run the dedicated binary

The verified local executable runs directly:

```sh
ISOGENY_PERF=/Volumes/SSD990-2/crypto/isogeny-performance-builds-20261014/rung3/release/isogeny-perf
"$ISOGENY_PERF" --help
"$ISOGENY_PERF" --out ./fresh-bmss-counts.jsonl \
  --revision f1321f6e6550b80bfafb44bd7ae1cd30816ab10c \
  --bits 192,224,256,511 --ells 7,31,101
"$ISOGENY_PERF" compare baseline.jsonl candidate.jsonl fresh-comparison.json
```

The executable includes its algorithm implementation and fixtures. It does not
read the repository or use Cargo at runtime. The revision argument records the
launcher-supplied source SHA; the executable hashes itself independently.
The file must be fresh. A missing PASS terminal record makes a run incomplete.

Retained baseline/candidate sources and binaries are under `evidence/baseline/`
and `evidence/rung1/`, `rung2/`, `rung3/`. The corrected baseline source is
`f7d5a9c198807c7ccf28979fd8c6ea6e60f10e35`, with binary SHA256
`9cb036ef3aa4998d722d4c43ed6d19ab33ee96c05ae98a5a3b18ed7c86abd245`.
`source-verified.tar.gz` is the complete baseline source archive; the initial
`source.tar.gz` attempt is retained as an interrupted archive.

`isogeny-perf campaign BASE CANDIDATE BASE_SHA CANDIDATE_SHA FRESH_DIR` runs
five A/A pairs and five interleaved A/B pairs at degrees 7/31. Append the native
`isolated-bench` path and reserved CPU list on Linux to collect isolation and
timing receipts. macOS collects repeated counts without wall samples.

The new ReleaseMe archive integration and four native runner profiles are in
`.github/workflows/release.yml` and `isogeny-performance.yml`. An actual asset
must be checked after publication; a successful local binary is not a
published release.

Maintainers build with Cargo. The native report helper links the same crate's
`json` and SHA256 implementation; it creates PDFs' Typst/SVG source from a PASS
comparison. Typst renders PDFs/150-dpi PNGs, and Poppler renders preview pages.
`make-report --manifest DIRECTORY FRESH_OUTPUT` creates a sorted byte/SHA256
inventory after all outputs are frozen. It rejects symlinks and special files.
