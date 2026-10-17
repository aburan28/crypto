#!/usr/bin/env python3
"""Lay out the finalized N83 evaluation report; never analyze timings.

Run only after the timed panel has stopped. The content must contain the native
summary's own values and intervals. This renderer does no statistical work.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def unresolved(value, path="content"):
    if isinstance(value, str):
        return [path] if "PENDING:" in value else []
    if isinstance(value, dict):
        return [p for key, item in value.items() for p in unresolved(item, f"{path}.{key}")]
    if isinstance(value, list):
        return [p for index, item in enumerate(value) for p in unresolved(item, f"{path}[{index}]")]
    return []


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--content", type=Path, default=ROOT / "report_content.json")
    parser.add_argument("--output", type=Path, default=ROOT / "n83_strong_reference_evaluation.pdf")
    parser.add_argument("--markdown", type=Path, default=ROOT / "report_source.md")
    parser.add_argument("--qa-dir", type=Path, default=ROOT / "preview")
    parser.add_argument("--allow-draft", action="store_true",
                        help="Permit explicitly marked pending facts in a draft")
    args = parser.parse_args()
    content = json.loads(args.content.read_text(encoding="utf-8"))
    missing = unresolved(content)
    if missing and not args.allow_draft:
        raise SystemExit("Native results or receipts are still pending at: " + ", ".join(missing))
    # Importing the engine, building and rasterizing occur only when invoked.
    # No benchmark values are recomputed or inferred here.
    from render_engine import build_pdf, quality_check, write_markdown
    metadata = build_pdf(args.content, args.output)
    write_markdown(metadata, args.markdown)
    qa = quality_check(args.output, args.qa_dir)
    print(json.dumps(qa, indent=2, ensure_ascii=False))
    if qa["page_bounds_violations"] or qa["replacement_glyphs"] or qa["missing_font_glyphs"]:
        raise SystemExit("Layout QA found an issue; inspect preview/layout_qa.json")


if __name__ == "__main__":
    main()
