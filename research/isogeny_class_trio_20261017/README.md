# Isogeny-class trio (2026-10-17)

Three roots for an experiment that iterates over every curve in an isogeny
class, one per standardized shape, each with Frobenius conductor `f > 1`:

| Shape | Curve | `f` | Strata | Class size |
|:--|:--|:--|:--|:--|
| Koblitz (K) | `icv1-f2m41-tm2308219-7f48b14a` | `409 · 1721` | 4 | 708,153 |
| Binary (B) | `icv1-f2m47-t4276051-826c5755` | `31^2` | 3 | 11,754,700 |
| Prime (P) | `icv1-fp47-tm8785349-43e36ce5` | `3^5` | 6 | 13,461,228 |

61-bit alternates: `icv1-f2m61-t158598901-ab42b6c5`,
`icv1-f2m61-tm2657393969-b9ff0fb0`, `icv1-fp61-tm813574577-a3b71255`.

- [Report](REPORT.md) and [PDF](REPORT.pdf); [frozen protocol](PROTOCOL.md)
- Figures: [strata](figures/strata.svg), [volcanoes](figures/volcanoes.svg),
  [search outcomes](figures/search-depth.svg), [Koblitz table](figures/koblitz-table.svg)
  (source [`figures/plot.py`](figures/plot.py))
- Records: [`data/`](data/), curve records [`curves.json`](curves.json)
- Native tools: [`selector/`](selector/) (search, class invariants, level tests;
  `cargo test --release`), [`identity/`](identity/) (ICV1 via `crypto_lib`)

PDF build: `pandoc -f gfm -t typst REPORT.md -o report.typ`, then compile
`report-main.typ` with Typst (the image lines are rewritten to
`#image(..., width: 100%)`).
