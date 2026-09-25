# Frozen row-certificate evidence

The preregistered source/input head is `ad2e447a15c4e57b2d8aabf7938646b83c1e9a2d`.
`FROZEN.json` SHA-256 is `e08fd485952cc8a3674964ed9b9d95f0f8fa1ba6978b9c0c8cc2713974eb5ae4`.
Its immutable #762 input archive SHA-256 is
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`.
The SHA-fixed planted input is a public synthetic test fixture, not a challenge target.

The accepted run produced `raw.tar.gz` (1,185,760 bytes), SHA-256
`863f50ea872ce01b1999d0ef67bc0341fa0c764f2ad4187d16736768b0ab9b3d`,
and `receipt.json`, SHA-256
`f88075ffbfb278cf9f584a3c3f139293856b27bca99c64dd964fe7995e3d0dcc`.
The receipt records UTC start/end, producer and verifier child costs, all raw
file hashes and the frozen manifest. It retains failure details and partial
file hashes if execution aborts. `input_construction_receipt.json` is frozen
with the source and records the bounded, byte-identical planted-input rerun.

`raw.tar.gz` contains the complete 2,003-target coset roster, all 17,917
positive n=13 rows, eight exceptional controls, eight public synthetic n=131
rows, per-cell summaries and `verify_report.json`. Each row contains source
factor indices/points, exact `Q,T`, transported/projected points, sign and
coefficient for every term, and the canonical aggregated row. A miss is
explicit in the roster. Independent replay reconstructs every row from the
original #762 archive and input fixture; it does not trust saved coefficients.

From this directory run `python3 ci_replay.py`. It checks the protocol-pinned
manifest, every frozen source/input hash, the older archive hash, recreates
the planted input, checks the output archive/receipt/raw hashes, then fully
replays the rows with the separate bit-serial/Fermat arithmetic from #762.
No producer execution is required for archive-only replay.
