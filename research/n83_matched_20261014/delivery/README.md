# Evidence packaging and delivery

`package_evidence.py` only packages existing files and checks their byte
identities. It does not calculate scientific statistics, run a solver, replay
arithmetic, select parameters or reinterpret a failure.

Run it after the native results, independent review, PDF inspection, benchmark
copies and canonical diagnostic pages are final. It refuses to overwrite an
existing final archive. `reports/run-20261014T1940Z/HASH_MANIFEST.json` records
each payload's byte count and SHA256. The archive includes that manifest, but
the manifest excludes itself and the external archive receipt to avoid a
circular identity. `EVIDENCE_ARCHIVE_RECEIPT.json` records the archive and
manifest identities, successful CRC checks and payload-hash checks. The ZIP
is copied unchanged to the committed report directory and `/workspace/report_output`.

The previous139-launch study is included as its unchanged historical evidence
ZIP; its own2854 entries receive a separate CRC check before packaging. Its
earlier source/fixture trees remain inside that archive, not rewritten as new
measurements. The original70-page comprehensive PDF is included separately
with its pinned SHA256. The new archive additionally contains this study's
raw campaign, verification attempts, native postprocessing sources and JSON,
presentation sources, inspected PDFs, graphs, benchmark copies, canonical
diagnostic snapshots and exact native executable files. Large build caches,
toolchain installations and regenerable page-preview PNGs are not included.
Graph PNG/SVG files, final PDFs, rendering sources and page-by-page QA receipts
are retained.

The campaign's native files retain original absolute paths and recorded
commands. Replaying those commands requires the recorded workspace layout
and compatible native dependencies; extracting an archive does not establish
an independent-host timing reproduction or a full compiled dependency closure.
The native checker validates answers, witnesses and file/matrix bindings;
the PDF and ZIP checks validate presentation and packaging, not mathematics.

Required human review and applicable current-head CI remain PR gates. Neither
a packaging PASS nor an inspected PDF authorizes bypassing those gates.
