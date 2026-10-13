#!/usr/bin/env python3
"""Render the SEC 2 inventory as a readable report, SVG and Typst source."""

from __future__ import annotations

import collections
import html
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
data = json.loads((HERE / "inventory.json").read_text())
rows = data["curves"]


def current(row: dict) -> bool:
    return "2010_v2" in row["listed_in"]


def field(row: dict) -> str:
    f = row["traits"]["field"]
    return ("Fp" if f["kind"] == "prime" else "F2m") + f"/{f['size_bits']}"


def cell(row: dict, key: str) -> str:
    t = row["traits"]
    if key == "cm":
        disc = t["frobenius"]["cm_disc"]
        if disc in ("-3", "-7", "-15"):
            return f"{disc} ({t['frobenius']['status']})"
        return t["frobenius"]["status"]
    if key == "subfield":
        s = t["subfield"]
        if t["field"]["kind"] == "prime":
            return "n/a"
        k = s["definition_degree"]
        return f"2^{k}" if k and k < t["field"]["degree"] else "full"
    if key == "embedding":
        e = t["embedding"]
        return "k>1000" if e["degree"] is None and e["lower_bound"] == 1000 else e["status"]
    if key == "cheon":
        a = row["cheon_n_pm_1_advisory"]
        return "/".join("F" if a[k]["status"] == "fails_advisory_bound" else "?"
                        for k in ("n_minus_1", "n_plus_1"))
    raise KeyError(key)


def render_svg() -> None:
    width, height = 950, 1090
    x0, x1 = 240, 825
    y0, dy = 92, 28
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="#ffffff"/>',
           '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#23313d}'
           '.title{font-size:22px;font-weight:700}.axis{font-size:12px;fill:#526573}'
           '.name{font-size:14px;font-weight:600}.small{font-size:11px;fill:#526573}</style>',
           '<text x="28" y="35" class="title">SEC 2 curve sizes and edition membership</text>',
           '<text x="28" y="57" class="small">Source: SEC 2 v1.0 (2000) and v2.0 (2010); field bits, 33 named curves</text>']
    for bits in (0, 112, 192, 256, 384, 521, 571):
        x = x0 + (x1 - x0) * bits / 571
        out += [f'<line x1="{x:.1f}" x2="{x:.1f}" y1="78" y2="1014" stroke="#e4e9ee"/>',
                f'<text x="{x:.1f}" y="1035" text-anchor="middle" class="axis">{bits}</text>']
    for i, row in enumerate(rows):
        y = y0 + i * dy
        bits = row["traits"]["field"]["size_bits"]
        x = x0 + (x1 - x0) * bits / 571
        color = "#176a9f" if current(row) else "#b86b16"
        out.append(f'<text x="28" y="{y+4}" class="name">{html.escape(row["name"])}</text>')
        out.append(f'<line x1="{x0}" x2="{x:.1f}" y1="{y}" y2="{y}" '
                   f'stroke="{color}" stroke-width="7" stroke-linecap="round"/>')
        if row["traits"]["field"]["kind"] == "prime":
            out.append(f'<circle cx="{x:.1f}" cy="{y}" r="6" fill="{color}"/>')
        else:
            out.append(f'<rect x="{x-6:.1f}" y="{y-6}" width="12" height="12" fill="{color}"/>')
        h = row["traits"]["subgroup"]["cofactor"]
        out.append(f'<text x="{x+12:.1f}" y="{y+4}" class="small">{bits} bits; h={h}</text>')
    out += ['<circle cx="32" cy="1060" r="6" fill="#176a9f"/>',
            '<text x="46" y="1064" class="small">2010 recommended; circle=prime field, square=binary field</text>',
            '<rect x="510" y="1054" width="12" height="12" fill="#b86b16"/>',
            '<text x="528" y="1064" class="small">2000 edition only</text>', '</svg>']
    (HERE / "sizes.svg").write_text("\n".join(out) + "\n")


def render_markdown() -> None:
    order_status = collections.Counter(r["traits"]["order_check"]["status"] for r in rows)
    cm_status = collections.Counter(r["traits"]["frobenius"]["status"] for r in rows)
    twist_status = collections.Counter(r["traits"]["twist"]["status"] for r in rows)
    fails = {k: sum(r["cheon_n_pm_1_advisory"][k]["status"] == "fails_advisory_bound" for r in rows)
             for k in ("n_minus_1", "n_plus_1")}
    out = ["# SEC 2 curve inventory and structural traits (2026-10-12)", "",
           "## Question and scope", "",
           "Inventory every named curve in SEC 2 v2.0 (2010), retain SEC 2 v1.0 (2000)"
           " entries omitted from v2.0, and bind each name to the `crypto` repository's exact"
           " curve identity and complete trait record. This is a catalog and arithmetic check,"
           " not an implementation or attack-cost measurement.", "",
           "## Sources and method", "",
           f"- **Primary:** [SEC 2 v1.0]({data['sources']['2000_v1']['url']})"
           f" (`{data['sources']['2000_v1']['sha256']}`) and"
           f" [SEC 2 v2.0]({data['sources']['2010_v2']['url']})"
           f" (`{data['sources']['2010_v2']['sha256']}`). Immutable PDF copies are in `sources/`.",
           f"- **Catalog:** `docs/curves/registry.json` SHA-256 `{data['catalog']['registry_sha256']}`;"
           f" `docs/curves/traits.json` SHA-256 `{data['catalog']['traits_sha256']}`."
           f" `docs/curves/covers.json` SHA-256 `{data['catalog']['covers_sha256']}`;"
           f" `docs/curves/jacobian-transfers.json` SHA-256"
           f" `{data['catalog']['transfers_sha256']}`. The source section pointers,"
           " curve identities, primary parameter checks, catalog traits, cover records,"
           " and transfer screening records are in `inventory.json`.",
           "- **Snapshot alignment:** `covers.json` is pinned to the current registry."
           " `traits.json` and `jacobian-transfers.json` carry older whole-registry"
           " hashes, but all 33 selected rows match the current exact slug and order/trace"
           " or ICV1, respectively. This report does not claim global freshness for"
           " those two derived files.",
           "- **Derived check:** `build.py` reads both PDFs through `pdftotext -layout`,"
           " compares field prime or reduction polynomial, coefficients, uncompressed base"
           " point, prime subgroup order, and cofactor with the exact registry representation,"
           " and refuses a mismatch. It extracts each stated seed without independently"
           " re-running the ANSI X9.62 seed-to-curve procedure. The existing `curve_traits`"
           " engine supplies CM, conductor, endomorphism bounds, subfields, embedding degree,"
           " twist, and small-isogeny local data, with each status preserved.",
           "- **Bounded SEC 1 check:** trial division by every prime up to 100,000 on"
           " `n−1` and `n+1`. If the remaining cofactor is at most `n^(19/20)`, no prime"
           " factor can satisfy the later [SEC 1 v2.0](https://www.secg.org/sec1-v2.pdf)"
           " advisory threshold. Other rows stay `undetermined`; the residual is not"
           " declared prime.", "",
           "## Results", "",
           "- **Primary match:** 33 distinct v1.0 names and 20 v2.0 names; all 53"
           " name-by-edition core parameter blocks agree with their catalog entry."
           " There are 15 prime-field and 18 binary-field names in v1.0;"
           " 8 prime-field and 12 binary-field names in v2.0.",
           "- **Historical entries:** " + ", ".join(f"`{r['name']}`" for r in rows if not current(r)) + ".",
           f"- **Trait status:** subgroup-order certificates {dict(order_status)};"
           f" CM/conductor {dict(cm_status)}; twist factorization {dict(twist_status)}."
           " Statuses are the catalog's and are not upgraded by this report.",
           "- **Covers and transfer screening:** all 33 curves have a matched"
           " `verified_over_declared_field` cover record and a matched Jacobian-transfer"
           " screening row. The cover records mark subgroup transport `not_tested` on"
           " all 33. The transfer rows all have verdict `boundary`: their formulas are"
           " screening models, not measured ECDLP gains.",
           f"- **SEC 1 n±1 advisory:** the bounded scan establishes failure of the"
           f" 19/20-factor recommendation for {fails['n_minus_1']} of 33 `n−1` values and"
           f" {fails['n_plus_1']} of 33 `n+1` values. The rest are undetermined under"
           " this trial limit. This is conditional Diffie–Hellman-oracle guidance,"
           " not a generic discrete-log cost result. The bounded check also finds"
           " 15 `n−1` and 10 `n+1` failures among the 20 curves retained in 2010;"
           " failure of this advisory condition alone cannot account for the 193-bit omission.", "",
           "![Field size by curve and edition](sizes.svg)", "",
           "### Curve-by-curve index", "",
           "The published base point has order `n`, and the whole curve group has"
           " order `#E(Fq) = h × n`. `h` is the cofactor. The table gives the bit length"
           " of `n`; its exact value and the generator coordinates are in `inventory.json`.", "",
           "`F/F` in the last column means the bounded scan established failure of"
           " the SEC 1 recommendation for both `n−1` and `n+1`; `?` means undetermined."
           " All exact integers, small-prime splitting, volcano depths, generator coordinates,"
           " seed bytes, and status-bearing trait objects are in `inventory.json`.", "",
           "| Curve | SEC 2 | Field | n bits | h | CM | End(E) | Subfield | Embedding | Twist | n−1/n+1 |",
           "| --- | --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- |"]
    for r in rows:
        t = r["traits"]
        out.append("| " + " | ".join([f"`{r['name']}`", "2010" if current(r) else "2000 only",
                                  field(r), str(t["subgroup"]["bits"]),
                                  t["subgroup"]["cofactor"], cell(r, "cm"),
                                  t["endomorphism"]["status"], cell(r, "subfield"),
                                  cell(r, "embedding"), t["twist"]["status"],
                                  cell(r, "cheon")]) + " |")
    out += ["", "## Interpretation and limits", "",
            "The standard's `r` and `k` name suffixes label its random-parameter and"
            " Koblitz families; they do not certify a security ranking. The seed bytes in"
            " `inventory.json` establish published provenance only. No curve-specific"
            " vulnerability, isogeny route, subgroup transport, or performance advantage"
            " follows from this inventory. Its `probable`, `bounded`, `unknown`, and"
            " `not_applicable` statuses remain distinct.", "",
            "The cover and transfer rows do not establish a subgroup-preserving"
            " usable route or a solver advantage. Their exact hypotheses and the"
            " absent transport checks remain visible in `inventory.json`.", "",
            "The current `crypto-autoresearcher:analysis/binstd-curve-audit/README.md`"
            " supplies an independent, secondary-source subfield scan for the binary"
            " names. Its scope is structural arithmetic and is kept separate from this"
            " primary SEC 2 parameter crosscheck.", "",
            "## Next check", "",
            "For an exact endomorphism ring or twist large-prime conclusion on a row"
            " still marked `bounded` or `unknown`, run a separately specified certificate"
            " or exact factorization for that row and attach it to its existing curve UID."
            " The present scan's unresolved residuals require either a certified prime"
            " factor exceeding `n^(19/20)` or enough additional factors to bound every"
            " residual below it. No experiment or research-state transition was made.", "",
            "## Reproduce", "",
            "From the repository root, run `python3 research/sec2_curve_inventory_20261012/build.py --check`."
            " To regenerate derived artifacts, run that script without `--check`, then"
            " `python3 research/sec2_curve_inventory_20261012/render.py`, then"
            " `typst compile research/sec2_curve_inventory_20261012/REPORT.typ"
            " research/sec2_curve_inventory_20261012/REPORT.pdf`. This rendering used"
            " Poppler `pdftotext` 26.04.0 and Typst 0.15.1.", ""]
    (HERE / "REPORT.md").write_text("\n".join(out))


def typst_escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("#", "\\#")


def render_typst() -> None:
    lines = [
        '#set page(paper: "a4", margin: 17mm)',
        '#set text(font: "Arial", size: 9pt)',
        '#set par(justify: true)',
        '= SEC 2 curves: identity and structural traits',
        '2026-10-12 · primary-standard crosscheck · 33 named curves',
        '',
        'The 2000 edition lists 33 curves. The 2010 edition retains 20 of them; '
        '13 are historical entries. Every field, coefficient, uncompressed base point, '
        'subgroup order, and cofactor in both editions matched the exact crypto registry record.',
        '',
        '#figure(image("sizes.svg", width: 100%), caption: [Field size and edition membership. '
        'Circle: prime field; square: binary field. Blue: listed in 2010. Amber: 2000 only.])',
        '',
        '#pagebreak()',
        '== Curve index',
        'Each published base point has order n, and the full curve has order h × n. '
        'Full exact identities, generators, n values, seeds, small-prime splitting, conductors, '
        'embedding and twist records, cover and transfer screening, and all status fields '
        'are in `inventory.json`. '
        'F in the last column means the SEC 1 advisory bound fails; ? means the bounded '
        'n±1 check was inconclusive.',
        '',
        '#set text(size: 7.3pt)',
        '#table(columns: (1.3fr, 0.8fr, 0.8fr, 0.5fr, 0.9fr, 0.8fr, 0.9fr), '
        'inset: 3pt, stroke: 0.35pt + rgb("dddddd"),',
        '[Curve], [Edition], [Field], [h], [CM status], [Twist], [n±1],',
    ]
    for r in rows:
        t = r["traits"]
        cells = [r["name"], "2010" if current(r) else "2000 only", field(r),
                 t["subgroup"]["cofactor"], cell(r, "cm"), t["twist"]["status"], cell(r, "cheon")]
        lines.append("  " + ", ".join(f'[{typst_escape(x)}]' for x in cells) + ",")
    lines += [')', '', '== Sources, method, and status',
              'SEC 2 v1.0 (2000): `https://www.secg.org/SEC2-Ver-1.0.pdf`',
              f"SHA-256 `{data['sources']['2000_v1']['sha256']}`.",
              'SEC 2 v2.0 (2010): `https://www.secg.org/sec2-v2.pdf`',
              f"SHA-256 `{data['sources']['2010_v2']['sha256']}`.", '',
              'The two PDFs supply 33 and 20 exact named parameter blocks. `build.py` compares '
              'field prime or reduction polynomial, coefficients, uncompressed base point, '
              'prime subgroup order, and cofactor with the registry. All 53 blocks match. '
              'Seeds are transcribed but the ANSI seed-to-curve derivation is not independently replayed.',
              '',
              'Existing `curve_traits` records give subgroup, CM/conductor, endomorphism bounds, '
              'small-prime splitting, subfield, embedding, and twist data. Cover and Jacobian-transfer '
              'screening rows match all 33 exact identities. Every cover has subgroup transport '
              '`not_tested`; the transfer formulas are screening models.', '',
              'The bounded n±1 scan certifies failure of the later SEC 1 19/20-factor advisory '
              'for 25 n−1 values and 19 n+1 values. Among curves retained in 2010, 15 n−1 '
              'and 10 n+1 values also fail, so failure of this advisory condition alone cannot account for the 193-bit omission. '
              'The other cases remain unresolved. '
              'The trait and transfer files have older whole-registry pins, but all 33 selected '
              'rows match current exact identities and order/trace.', '',
              '== Evidence boundaries',
              'This is a primary-parameter check and a structured view of existing catalog traits. '
              'A `probable` or `bounded` trait is not promoted to proved. The n±1 scan uses trial '
              'division through 100,000 and proves failure of the SEC 1 advisory threshold only '
              'when the residual is already below it. Covers and transfer screening do not verify '
              'subgroup transport or measure attack cost. Full status-bearing records are in '
              '`inventory.json`; source PDF copies are in `sources/`.']
    (HERE / "REPORT.typ").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    render_svg()
    render_markdown()
    render_typst()
    print("Rendered REPORT.md, REPORT.typ, and sizes.svg")
