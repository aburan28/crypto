"""Independent PDF presentation QA only; no research-result computation."""
from pathlib import Path
import hashlib
import json
import fitz

ROOT = Path(__file__).resolve().parent
RUN = ROOT.parent
ADDENDUM = RUN / "n83_matched_ic_rho_evaluation.pdf"
COMBINED = RUN / "unified_curve_basis_report_matched_ic_rho.pdf"
PRIOR = Path("/workspace/report_output/unified_curve_basis_report_faster_rho.pdf")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def inspect(path):
    doc = fitz.open(path)
    bounds, glyphs, links, invalid_destinations = [], [], [], []
    texts = []
    for index, page in enumerate(doc):
        texts.append(page.get_text())
        for link in page.get_links():
            item = {"page": index + 1, "kind": link["kind"],
                    "uri": link.get("uri"), "destination_page": link.get("page")}
            links.append(item)
            if link["kind"] == fitz.LINK_GOTO and not 0 <= link.get("page", -1) < len(doc):
                invalid_destinations.append(item)
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line["spans"]:
                    box = fitz.Rect(span["bbox"])
                    if (box.x0 < -0.5 or box.y0 < -0.5 or
                            box.x1 > page.rect.width + 0.5 or
                            box.y1 > page.rect.height + 0.5):
                        bounds.append({"page": index + 1, "bbox": list(box),
                                       "text": span["text"]})
                    if "\ufffd" in span["text"] or "\x00" in span["text"]:
                        glyphs.append({"page": index + 1, "text": span["text"]})
    return doc, {"file": str(path), "sha256": digest(path.read_bytes()),
                 "pages": len(doc), "links": links,
                 "link_annotations": len(links), "bounds_violations": bounds,
                 "replacement_or_null_glyphs": glyphs,
                 "invalid_internal_link_destinations": invalid_destinations}, texts


def compare_pages(source, target, offset, scale=2):
    results = []
    for index, page in enumerate(source):
        appended = target[index + offset]
        left = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
        right = appended.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False)
        results.append({"source_page": index + 1,
                        "combined_page": index + offset + 1,
                        "text_equal": page.get_text() == appended.get_text(),
                        "size_equal": (left.width, left.height, left.n) ==
                                      (right.width, right.height, right.n),
                        "dpi": 72 * scale,
                        "pixels_equal": left.samples == right.samples,
                        "source_pixels_sha256": digest(left.samples),
                        "combined_pixels_sha256": digest(right.samples)})
    return results


add, add_info, add_text = inspect(ADDENDUM)
combined, combined_info, _ = inspect(COMBINED)
prior, prior_info, _ = inspect(PRIOR)
assert len(add) == 17 and len(prior) == 70 and len(combined) == 87
historical = compare_pages(prior, combined, len(add))
new_pages = compare_pages(add, combined, 0)
new_pages_72 = compare_pages(add, combined, 0, scale=1)
all_preserved = lambda rows: all(r["text_equal"] and r["size_equal"] and
                                 r["pixels_equal"] for r in rows)
for index, page in enumerate(add):
    page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False).save(
        str(ROOT / f"page_{index + 1:02d}.png"))
text = "\n\n".join(f"PAGE {index + 1}\n{value}" for index, value in enumerate(add_text))
(ROOT / "addendum_text.txt").write_text(text)
visual_path = ROOT / "visual_review.json"
visual = json.loads(visual_path.read_text()) if visual_path.exists() else {}
visual_bound = (visual.get("status") == "PASS" and
                visual.get("addendum_sha256") == add_info["sha256"] and
                visual.get("combined_sha256") == combined_info["sha256"])
receipt = {"kind": "independent_pdf_layout_and_preservation_qa",
           "research_results_recomputed": False,
           "conductor_check": {"status": "unavailable", "exit_code": 127,
                               "reason": "conductor: command not found"},
           "addendum": add_info, "combined": combined_info, "historical": prior_info,
           "page_count_sum_matches": len(combined) == len(add) + len(prior),
           "historical_pages_at_144dpi": historical,
           "new_pages_at_144dpi": new_pages,
           "new_pages_at_72dpi": new_pages_72,
           "all_70_historical_pages_preserved": all_preserved(historical),
           "all_17_new_pages_preserved_in_combined": all_preserved(new_pages),
           "all_17_new_pages_preserved_in_combined_at_72dpi": all_preserved(new_pages_72),
           "pending_marker_present": "PENDING:" in text,
           "visual_review_status": "PASS" if visual_bound else "pending_independent_page_by_page_inspection",
           "visual_review_receipt": str(visual_path) if visual_bound else None,
           "external_link_scope": "Annotation destinations inspected; remote HTTP availability not fetched."}
(ROOT / "qa.json").write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps({"pages": [len(add), len(prior), len(combined)],
                  "new_link_annotations": add_info["link_annotations"],
                  "new_bounds_violations": len(add_info["bounds_violations"]),
                  "new_glyph_issues": len(add_info["replacement_or_null_glyphs"]),
                  "invalid_internal_destinations": len(add_info["invalid_internal_link_destinations"]),
                  "prior70_text_and_pixels_144dpi": all_preserved(historical),
                  "new17_text_and_pixels_144dpi": all_preserved(new_pages),
                  "new17_text_and_pixels_72dpi": all_preserved(new_pages_72),
                  "pending_marker_present": receipt["pending_marker_present"]}, indent=2))
