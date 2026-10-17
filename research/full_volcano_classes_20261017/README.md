# Three additional complete volcano graphs

The three additional classes contain 24,823 vertices, 69,300 directed
prime-degree maps and 51 conductor strata. Native construction, full map
replay, copied-database replay and exact stratum populations passed. These
are tractable additions; they do not complete the original A/C/D/E runs.

| Source model | Prime-field bits | Maximum End(E) conductor | Vertices | Maps | Strata | Degrees |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| icv1-fp57-t536870914-6cabb0c9 | 57 | 8192 | 8192 | 16384 | 14 | 2 |
| icv1-fp57-tm536870941-e6d0c766 | 57 | 6561 | 6561 | 13122 | 9 | 3 |
| icv1-fp59-t1073741954-4dd0b1ef | 59 | 1728 | 10070 | 39794 | 28 | 2,3 |

Source inputs begin at conductor 1. The `inputs/*-floor.json` models have
the listed maximum End(E) conductors, established by the complete labeled
graph rather than equating the source conductor with the Frobenius conductor.
The archive contains SQLite databases, all-node CSVs, all-kernel JSONL,
stratum adjacency counts and an explicit-map example for every stratum.

[Download the complete archive](https://drive.google.com/file/d/1SSyQ5bY0RPTK_abRW6MCG78cdANonBVi/view?usp=drivesdk)
or [open the PDF](https://drive.google.com/file/d/1E0QPO5jRkOJ9iFD4EQvWymeguxp_QVyI/view?usp=drivesdk).
Drive metadata readback verified both uploads' names, sizes and parent folder.
The local ZIP CRC and every delivered file's SHA-256 passed; a remote whole-file
hash was not checked. The 10,066,274-byte archive SHA-256 is
`17c2d5d3290925b5fe0d84534423dd6e85b1804d4382dad0c00955b47a76458e`.

See [the report](reports/run-20261017/REPORT.md), [PDF](reports/run-20261017/report.pdf),
[native counts](evidence/coverage.json), [protocol](PROTOCOL.md) and
`delivery.json` for hashes and verified Drive links.

```sh
# From the downloadable archive:
./bin/isogeny-class status --input graphs/mixed.sqlite
./bin/isogeny-class vertex --input graphs/mixed.sqlite --id 0
# verify refreshes metadata: replay a copy of the byte-frozen original.
cp graphs/mixed.sqlite mixed-replay.sqlite
./bin/isogeny-class verify --input mixed-replay.sqlite
# Generate the complete class starting at the largest-conductor model:
./bin/isogeny-class build --curve inputs/mixed-floor.json --out fresh.sqlite --steps all
```

Generation used protocol commit `3b7f4b39c` and frozen input/run commit
`946ab7419`. Source fingerprint and delivered binary hash are retained in
`evidence/implementation.sha256`. The executable/source snapshot is the
previously delivered `7c2e8f46a` implementation. New native fixture and export
sources are included separately. Run the repository version with
`bash run.sh /path/to/isogeny-class /fresh/graph/directory two`, then repeat
for `three` and `mixed`; output paths must not already contain a graph.

The 6 native unit and 2 process integration tests passed again on the updated
checkout. Full replay shares arithmetic with construction; independent
implementation replay remains NOT_RUN. Graph-wide ECDLP and index-calculus
timing, operation counts, min/max/mean and relation yield remain NOT_RUN.
All-degree isogenies form an infinite family; only the displayed finite sets
of prime degrees are enumerated here.

The first metadata import failed because stdout was redirected to an empty
`.json` file inside the directory being scanned. The empty file and error
were retained; using a `.txt` stdout log fixed the orchestration. No graph
data or cryptographic arithmetic changed. Conductor reservation was attempted
but unavailable (127). Canonical IC/rho timing graphs are unchanged because
this run supplies graph coverage, not solver performance measurements.

The catalogue cover generator and full recheck verified all 2,126 registered
models. Browser and leaderboard generation/freshness checks passed. The full
debug traits rebuild asserted in existing binary-subfield code and was stopped;
its queued full check was not run. A native incremental refresh computed the
six new prime-model records, preserved all 2,296 baseline rows exactly, and
verified unchanged hashes for all other inputs. This is explicitly not a full
traits replay: see [the open validation problem](../../problems/curve-traits-debug-subfield-trace.md).
Global naming lint still reports five legacy occurrences in two unchanged
archive reports; frozen evidence was preserved.
