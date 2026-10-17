# Isogeny-class trio (2026-10-17)

Three roots for an experiment that iterates over every curve in an isogeny
class, one per standardized shape, each with Frobenius conductor `f > 1`.
The user chose the 61-bit trio:

| Shape | Curve | `f` | Strata | Class size | Root |
|:--|:--|:--|:--|:--|:--|
| Koblitz (K) | `icv1-f2m61-t158598901-ab42b6c5` | `1951 · 587551` | 4 | 1,147,491,009 | crater |
| Binary (B) | `icv1-f2m61-tm2657393969-b9ff0fb0` | `7^2` | 3 | 851,960,648 | floor |
| Prime (P), `p = 2^61 − 1` | `icv1-fp61-tm813574577-a3b71255` | `3^2 · 5^3` | 12 | 920,783,750 | floor |

The 41/47-bit trio (`icv1-f2m41-tm2308219-7f48b14a`,
`icv1-f2m47-t4276051-826c5755`, `icv1-fp47-tm8785349-43e36ce5`) stays in the
report for reference. All six are listed in the curve testbed,
[`docs/curves/testbed/`](../../docs/curves/testbed/README.md).

- [Report](REPORT.md) and [PDF](REPORT.pdf); [frozen protocol](PROTOCOL.md)
  with its amendment
- Figures: [strata](figures/strata.svg), [volcanoes](figures/volcanoes.svg),
  [search outcomes](figures/search-depth.svg), [Koblitz table](figures/koblitz-table.svg)
  (source [`figures/plot.py`](figures/plot.py))
- Records: [`data/`](data/), curve records [`curves.json`](curves.json)
- Native tools: [`selector/`](selector/) (search, class invariants, level tests;
  `cargo test --release`), [`identity/`](identity/) (ICV1 via `crypto_lib`)

PDF build: `pandoc -f gfm -t typst REPORT.md -o report.typ`, then compile
`report-main.typ` with Typst (the image lines are rewritten to
`#image(..., width: 100%)`).
