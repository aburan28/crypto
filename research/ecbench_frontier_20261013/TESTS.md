# Validation

Repository commit: `759978250c1ce53d7b87bd8a1c2304781fcc4fc4`.

| Check | Result |
|---|---|
| `cargo build --release --bin ecbench` under Rust 1.94.1 | pass; pre-existing warnings only |
| factor-base K17 native verify, 8 replays | pass, 8/8 identical; receipt `9306d33e…` |
| F4/F5/F6 K13 native verify, 7 replays | pass, 7/7 identical across all six arms; receipt `0304c496…` |
| root-finder K13 native verify, sampled 8 replays | pass, 8/8 identical; sampler covered rho rows; receipt `4814136e…` |
| root-finder K13 native verify, replay-all | **fail retained**, 120/168 identical; all 48 inherited-F4 rows differ in solver phase counts; receipt `6110727a…` |
| `jq empty frontier-data.json` | pass |
| `python -m py_compile .../plot.py` | pass |
| render both SVG and 150 dpi PNG figures from `frontier-data.json` | pass |
| Typst compile `REPORT.typ` | pass, 5-page A4 PDF, SHA-256 `4fdc1540fa20e66bb63da7f83cd00499f4c0f08e8f1d831512b286bc9c88a908` |
| visual inspection of all PDF pages | pass; tables, captions, identifiers, and graph labels legible |
| benchmark copies byte-equal to study PDF/SVG/PNG | pass by SHA-256 |

The failed replay is an evidence result, not a waived test. It is described in
the README and PDF and blocks promotion of the affected F4 stage counts.
