#!/usr/bin/env python3
"""Render a native-result addendum and prepend it to the frozen prior PDF.

Document rendering/merging only. No experiment, selector, statistics, interval
estimation, cryptographic arithmetic or benchmark result is computed here.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent


def pending(value, location="content"):
    if isinstance(value, str):
        return [location] if "PENDING:" in value else []
    if isinstance(value, dict):
        return [p for k, v in value.items() for p in pending(v, f"{location}.{k}")]
    if isinstance(value, list):
        return [p for n, v in enumerate(value) for p in pending(v, f"{location}[{n}]")]
    return []


def diagram_previews(content):
    for section in content.get("sections", []):
        figures = [section["figure"]] if section.get("figure") else []
        figures += [b for b in section.get("blocks", []) if b.get("type") == "figure"]
        for figure in figures:
            png = ROOT / figure["path"]
            svg = png.with_suffix(".svg")
            if svg.exists() and png.suffix.lower() == ".png":
                subprocess.run(["inkscape", str(svg), "--export-type=png",
                                f"--export-filename={png}", "--export-width=2200"], check=True)


def combine(addendum, prior, output, expected_prior_pages):
    import fitz
    addendum, prior, output = [Path(p).resolve() for p in (addendum, prior, output)]
    if output in (addendum, prior):
        raise ValueError("Combined output must be a new file, distinct from both source PDFs")
    with fitz.open(addendum) as first, fitz.open(prior) as historical:
        if len(historical) != expected_prior_pages:
            raise ValueError(f"Expected frozen prior PDF with {expected_prior_pages} pages; found {len(historical)}")
        result = fitz.open()
        result.insert_pdf(first)
        result.insert_pdf(historical)
        toc = first.get_toc()
        toc.append([1, "Historically frozen prior report — default-policy observations", len(first) + 1])
        toc += [[level + 1, title, page + len(first)] for level, title, page in historical.get_toc()]
        result.set_toc(toc)
        result.set_metadata({"title":"Unified curve-basis report with optimized-rho addendum",
                             "subject":"New tuning addendum followed by historically frozen prior report"})
        # This is page-preservation layout QA, not research-result analysis.
        preserved = [historical[n].get_text() == result[len(first)+n].get_text()
                     for n in range(len(historical))]
        if not all(preserved):
            raise ValueError("Historical page text was not preserved by the PDF merge")
        output.parent.mkdir(parents=True, exist_ok=True)
        result.save(output, garbage=4, deflate=True)
        receipt = {"kind":"pdf_layout_merge", "addendum":str(addendum),
                   "historical_source_read_only":str(prior),
                   "output":str(output), "addendum_pages":len(first),
                   "historical_pages":len(historical), "combined_pages":len(result),
                   "historical_page_text_preserved":True,
                   "research_results_recomputed":False}
        result.close()
        return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--content", type=Path, default=ROOT/"report_content.json")
    parser.add_argument("--output", type=Path, default=ROOT/"n83_optimized_rho_evaluation.pdf")
    parser.add_argument("--prior-report", type=Path,
                        default=Path("/workspace/report_output/unified_curve_basis_report_updated_n83.pdf"))
    parser.add_argument("--expected-prior-pages", type=int, default=48)
    parser.add_argument("--user-addendum", type=Path,
                        default=Path("/workspace/report_output/n83_optimized_rho_evaluation.pdf"))
    parser.add_argument("--combined-output", type=Path,
                        default=Path("/workspace/report_output/unified_curve_basis_report_optimized_rho.pdf"))
    parser.add_argument("--addendum-only", action="store_true")
    args = parser.parse_args()
    content = json.loads(args.content.read_text())
    missing = pending(content)
    if missing:
        raise SystemExit("Native findings are pending; render is withheld at: " + ", ".join(missing))
    # Deferred imports and figure jobs keep source preparation light while a
    # timed native campaign is active. Invoke only after root ends timing.
    from render_engine import build_pdf, quality_check, write_markdown
    diagram_previews(content)
    metadata = build_pdf(args.content, args.output)
    write_markdown(metadata, ROOT/"report_source.md")
    qa = quality_check(args.output, ROOT/"preview")
    if qa["page_bounds_violations"] or qa["replacement_glyphs"] or qa["missing_font_glyphs"]:
        raise SystemExit("Layout QA found an issue; inspect preview/layout_qa.json")
    args.user_addendum.parent.mkdir(parents=True, exist_ok=True)
    if args.output.resolve() != args.user_addendum.resolve():
        shutil.copy2(args.output, args.user_addendum)
    result = {"addendum_layout":qa}
    if not args.addendum_only:
        result["merge"] = combine(args.output,args.prior_report,args.combined_output,
                                  args.expected_prior_pages)
        result["combined_layout"] = quality_check(args.combined_output,ROOT/"combined_preview")
        if result["combined_layout"]["page_bounds_violations"] or result["combined_layout"]["replacement_glyphs"]:
            raise SystemExit("Combined document layout QA found an issue")
    (ROOT/"render_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
