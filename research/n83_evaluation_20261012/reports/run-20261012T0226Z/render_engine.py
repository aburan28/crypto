#!/usr/bin/env python3
"""Render the evidence synthesis JSON as a typeset, linked PDF.

This module only lays out supplied narrative, tables and figures. It does not
compute mathematical or benchmark results. Run with --help for paths/options.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import re
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Flowable, Frame, Image, KeepTogether, PageBreak,
    PageTemplate, Paragraph, Spacer, Table, TableStyle,
)
from reportlab.platypus.tableofcontents import TableOfContents


ROOT = Path(__file__).resolve().parent
PAGE_W, PAGE_H = A4
MARGIN = 45
CONTENT_W = PAGE_W - 2 * MARGIN
NAVY = colors.HexColor("#142B45")
TEAL = colors.HexColor("#087F83")
INK = colors.HexColor("#25354A")
MUTED = colors.HexColor("#5B6A7A")
RULE = colors.HexColor("#DAE4EB")
PALE = colors.HexColor("#EFF6F7")
MISSING_GLYPHS: set[str] = set()


def register_fonts():
    directory = Path("/usr/share/fonts/truetype/dejavu")
    names = {
        "DejaVuSans": "DejaVuSans.ttf",
        "DejaVuSans-Bold": "DejaVuSans-Bold.ttf",
        "DejaVuSans-Oblique": "DejaVuSans-Oblique.ttf",
        "DejaVuSans-BoldOblique": "DejaVuSans-BoldOblique.ttf",
        "DejaVuSerif": "DejaVuSerif.ttf",
        "DejaVuSerif-Bold": "DejaVuSerif-Bold.ttf",
        "DejaVuSerif-Italic": "DejaVuSerif-Italic.ttf",
        "DejaVuSerif-BoldItalic": "DejaVuSerif-BoldItalic.ttf",
        "DejaVuSansMono": "DejaVuSansMono.ttf",
        "DejaVuSansMono-Bold": "DejaVuSansMono-Bold.ttf",
        "DejaVuMath": "DejaVuMathTeXGyre.ttf",
    }
    for name, filename in names.items():
        pdfmetrics.registerFont(TTFont(name, str(directory / filename)))
    pdfmetrics.registerFont(TTFont("NotoPunctuation", "/usr/share/fonts/truetype/noto/NotoSansYi-Regular.ttf"))
    for family, italic in (("DejaVuSans", "Oblique"), ("DejaVuSerif", "Italic")):
        pdfmetrics.registerFontFamily(
            family, normal=family, bold=family + "-Bold",
            italic=family + "-" + italic, boldItalic=family + "-Bold" + italic,
        )
    pdfmetrics.registerFontFamily("DejaVuSansMono", normal="DejaVuSansMono",
                                  bold="DejaVuSansMono-Bold",
                                  italic="DejaVuSansMono",
                                  boldItalic="DejaVuSansMono-Bold")
    pdfmetrics.registerFontFamily("DejaVuMath", normal="DejaVuMath",
                                  bold="DejaVuMath", italic="DejaVuMath",
                                  boldItalic="DejaVuMath")
    pdfmetrics.registerFontFamily("NotoPunctuation", normal="NotoPunctuation",
                                  bold="NotoPunctuation", italic="NotoPunctuation",
                                  boldItalic="NotoPunctuation")


def make_styles():
    base = ParagraphStyle("Base", fontName="DejaVuSans", textColor=INK,
                          alignment=TA_LEFT)
    common = dict(parent=base)
    return {
        "body": ParagraphStyle("Body", **common, fontSize=10, leading=14.5,
                               spaceAfter=9, splitLongWords=1, allowWidows=0,
                               allowOrphans=0, uriWasteReduce=0.35),
        "section": ParagraphStyle("Section", **common, fontName="DejaVuSans-Bold",
                                  fontSize=21, leading=27, textColor=NAVY,
                                  spaceAfter=18, keepWithNext=1),
        "eyebrow": ParagraphStyle("Eyebrow", **common, fontName="DejaVuSans-Bold",
                                  fontSize=8.7, leading=13, textColor=TEAL,
                                  spaceAfter=9, tracking=1.1),
        "cover_title": ParagraphStyle("CoverTitle", fontName="DejaVuSerif",
                                      textColor=NAVY, fontSize=30, leading=38,
                                      spaceAfter=22),
        "cover_subtitle": ParagraphStyle("CoverSubtitle", **common,
                                         fontSize=13, leading=20,
                                         textColor=MUTED, spaceAfter=22),
        "cover_meta": ParagraphStyle("CoverMeta", **common, fontSize=9.5,
                                     leading=15, textColor=MUTED, spaceAfter=10),
        "callout": ParagraphStyle("Callout", **common, fontSize=10.3,
                                  leading=15.2, spaceAfter=0),
        "table": ParagraphStyle("TableCell", **common, fontSize=8.8,
                                leading=12, splitLongWords=1, spaceAfter=0,
                                uriWasteReduce=0.25),
        "table_head": ParagraphStyle("TableHead", fontName="DejaVuSans-Bold",
                                     textColor=colors.white, fontSize=8.8,
                                     leading=12, splitLongWords=1),
        "caption": ParagraphStyle("Caption", **common, fontSize=8.5,
                                  leading=12.2, textColor=MUTED, spaceAfter=12),
        "bullet": ParagraphStyle("Bullet", **common, fontSize=10, leading=14.5,
                                 leftIndent=13, firstLineIndent=-10,
                                 spaceAfter=7, splitLongWords=1),
        "small": ParagraphStyle("Small", **common, fontSize=8.3, leading=12,
                                spaceAfter=7, textColor=MUTED),
    }


_TAG = re.compile(r"(</?(?:b|i|u|sub|super|sup|font|link|a|br|strong|em|strike)(?:\s[^<>]*?)?\s*/?>)", re.I)


def rich_text(value):
    """Escape plain prose while preserving the supported Paragraph XML tags.

    Math symbols unsupported by the body font use DejaVu Math. Raw comparison
    signs and ampersands are safe. Plain newlines become line breaks.
    """
    body_chars = pdfmetrics.getFont("DejaVuSans").face.charToGlyph
    math_chars = pdfmetrics.getFont("DejaVuMath").face.charToGlyph
    punctuation_chars = pdfmetrics.getFont("NotoPunctuation").face.charToGlyph
    parts = _TAG.split(str(value))
    out = []
    for part in parts:
        if _TAG.fullmatch(part or ""):
            out.append(part.replace("<sup>", "<super>").replace("</sup>", "</super>"))
            continue
        # Decode once so already escaped authors' text is not double escaped.
        for character in html.unescape(part):
            if character == "\n":
                out.append("<br/>")
            elif character == "\t":
                out.append(" ")
            elif ord(character) in body_chars or character.isspace():
                out.append(html.escape(character, quote=False))
            elif ord(character) in math_chars:
                out.append('<font name="DejaVuMath">' + html.escape(character, quote=False) + "</font>")
            elif ord(character) in punctuation_chars:
                out.append('<font name="NotoPunctuation">' + html.escape(character, quote=False) + "</font>")
            else:
                MISSING_GLYPHS.add(character)
                out.append(html.escape(character, quote=False))
    return "".join(out)


def plain_text(value):
    return html.unescape(_TAG.sub("", str(value))).replace("\n", " ")


class AccentRule(Flowable):
    def __init__(self, width=CONTENT_W, height=12):
        super().__init__()
        self.width, self.height = width, height

    def draw(self):
        self.canv.setFillColor(TEAL)
        self.canv.rect(0, self.height - 3, 42, 3, fill=1, stroke=0)
        self.canv.setStrokeColor(RULE)
        self.canv.line(52, self.height - 1.5, self.width, self.height - 1.5)


class ReportDocument(BaseDocTemplate):
    def __init__(self, filename, metadata, **kwargs):
        super().__init__(filename, pagesize=A4, leftMargin=MARGIN,
                         rightMargin=MARGIN, topMargin=62, bottomMargin=49,
                         title=plain_text(metadata["title"]),
                         author=metadata.get("author", "Research synthesis"),
                         subject=plain_text(metadata.get("subtitle", "")),
                         pageCompression=1, **kwargs)
        self.metadata = metadata
        self.current_section = "Contents"
        self.section_page = None
        frame = Frame(MARGIN, self.bottomMargin, CONTENT_W,
                      PAGE_H - self.topMargin - self.bottomMargin,
                      leftPadding=0, rightPadding=0, topPadding=0,
                      bottomPadding=0, id="body")
        self.addPageTemplates(PageTemplate(id="report", frames=[frame],
                                           onPage=self.page_art))

    def page_art(self, canvas, doc):
        canvas.saveState()
        if doc.page == 1:
            canvas.setFillColor(NAVY)
            canvas.rect(0, PAGE_H - 17, PAGE_W, 17, fill=1, stroke=0)
            canvas.setFillColor(PALE)
            canvas.rect(PAGE_W - 150, 0, 150, 160, fill=1, stroke=0)
            canvas.setFillColor(TEAL)
            canvas.rect(PAGE_W - 150, 145, 70, 4, fill=1, stroke=0)
        else:
            canvas.setFillColor(MUTED)
            canvas.setFont("DejaVuSans", 7.5)
            header = self.metadata.get("running_header", "CURVE BASES  /  EVIDENCE SYNTHESIS")
            canvas.drawString(MARGIN, PAGE_H - 32, plain_text(header)[:92])
            canvas.setStrokeColor(RULE)
            canvas.setLineWidth(0.6)
            canvas.line(MARGIN, PAGE_H - 43, PAGE_W - MARGIN, PAGE_H - 43)
            canvas.line(MARGIN, 34, PAGE_W - MARGIN, 34)
            canvas.setFont("DejaVuSans", 7.3)
            canvas.drawString(MARGIN, 21, str(self.metadata.get("date", "")))
            canvas.drawRightString(PAGE_W - MARGIN, 21, f"{doc.page:02d}")
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if hasattr(flowable, "section_key"):
            title = flowable.getPlainText()
            key = flowable.section_key
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(title, key, level=0, closed=False)
            self.notify("TOCEntry", (0, title, self.page, key))


def render_table(spec, styles, default_padding=8):
    headers = spec.get("headers", [])
    rows = spec.get("rows", [])
    if not headers:
        raise ValueError("Each table needs a nonempty headers list")
    if any(len(row) != len(headers) for row in rows):
        raise ValueError("Every table row must have the same length as headers")
    widths = spec.get("widths", [1] * len(headers))
    if len(widths) != len(headers) or any(float(w) <= 0 for w in widths):
        raise ValueError("Table widths must be positive relative weights, one per column")
    total = sum(float(w) for w in widths)
    widths = [CONTENT_W * float(w) / total for w in widths]
    cells = [[Paragraph(rich_text(v), styles["table_head"]) for v in headers]]
    cells += [[Paragraph(rich_text(v), styles["table"]) for v in row] for row in rows]
    table = Table(cells, colWidths=widths, repeatRows=1, hAlign="LEFT",
                  splitByRow=1, splitInRow=1, spaceBefore=5, spaceAfter=15)
    padding = float(spec.get("padding", default_padding))
    commands = [
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), padding),
        ("BOTTOMPADDING", (0, 0), (-1, -1), padding),
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, TEAL),
        ("LINEBELOW", (0, 1), (-1, -1), 0.35, RULE),
    ]
    for i in range(1, len(cells)):
        if i % 2:
            commands.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#F4F7F9")))
    table.setStyle(TableStyle(commands))
    return table


def render_figure(spec, styles, source_dir):
    path = Path(spec["path"])
    if not path.is_absolute():
        path = source_dir / path
    if not path.exists():
        raise FileNotFoundError(f"Figure does not exist: {path}")
    picture = Image(str(path))
    max_height = float(spec.get("max_height", 365))
    factor = min(CONTENT_W / picture.imageWidth, max_height / picture.imageHeight)
    picture.drawWidth = picture.imageWidth * factor
    picture.drawHeight = picture.imageHeight * factor
    picture.hAlign = "CENTER"
    items = [Spacer(1, 6), picture, Spacer(1, 7)]
    if spec.get("caption"):
        items.append(Paragraph(rich_text(spec["caption"]), styles["caption"]))
    return KeepTogether(items)


def render_callout(value, styles):
    cell = Paragraph(rich_text(value), styles["callout"])
    box = Table([[cell]], colWidths=[CONTENT_W], spaceBefore=2, spaceAfter=15)
    box.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("LINEBEFORE", (0, 0), (0, -1), 3, TEAL),
        ("LEFTPADDING", (0, 0), (-1, -1), 13),
        ("RIGHTPADDING", (0, 0), (-1, -1), 13),
        ("TOPPADDING", (0, 0), (-1, -1), 12),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 12),
    ]))
    return box


def render_block(block, styles, source_dir):
    kind = block.get("type", "paragraph")
    if kind in ("paragraph", "text"):
        return [Paragraph(rich_text(block["text"]), styles["body"])]
    if kind == "small":
        return [Paragraph(rich_text(block["text"]), styles["small"])]
    if kind == "callout":
        return [render_callout(block["text"], styles)]
    if kind == "table":
        return [render_table(block, styles)]
    if kind == "figure":
        return [render_figure(block, styles, source_dir)]
    if kind == "bullets":
        return [Paragraph("• " + rich_text(v), styles["bullet"]) for v in block["items"]]
    if kind == "subheading":
        sub = ParagraphStyle("Subheading", parent=styles["body"],
                             fontName="DejaVuSans-Bold", fontSize=12.1,
                             leading=17, textColor=NAVY, spaceBefore=10,
                             spaceAfter=8, keepWithNext=1)
        return [Paragraph(rich_text(block["text"]), sub)]
    if kind == "page_break":
        return [PageBreak()]
    raise ValueError(f"Unknown block type: {kind}")


def build_pdf(content_path, output_path):
    register_fonts()
    content_path = Path(content_path).resolve()
    metadata = json.loads(content_path.read_text(encoding="utf-8"))
    styles = make_styles()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc = ReportDocument(str(output_path), metadata)
    story = [Spacer(1, 65),
             Paragraph(rich_text(metadata.get("cover_label", "RESEARCH SYNTHESIS")), styles["eyebrow"]),
             AccentRule(), Spacer(1, 20),
             Paragraph(rich_text(metadata["title"]), styles["cover_title"])]
    if metadata.get("subtitle"):
        story.append(Paragraph(rich_text(metadata["subtitle"]), styles["cover_subtitle"]))
    story.append(Spacer(1, 25))
    cover_meta = metadata.get("cover_meta", [str(metadata.get("date", "")),
                                             "Evidence, algorithms, and measurement protocol"])
    story.extend(Paragraph(rich_text(v), styles["cover_meta"]) for v in cover_meta if v)
    if metadata.get("cover_note"):
        story.extend([Spacer(1, 20), render_callout(metadata["cover_note"], styles)])
    if metadata.get("include_contents", True):
        story.extend([PageBreak(), Paragraph("Contents", styles["section"]),
                      Paragraph(rich_text(metadata.get("contents_intro", "A source-bound comparison of the exact fixture, recorded performance, correctness receipts, and admission limits.")), styles["small"]),
                      Spacer(1, 10)])
        toc = TableOfContents()
        toc.levelStyles = [ParagraphStyle(
            "ContentsEntry", fontName="DejaVuSans", fontSize=10.2,
            leading=15.5, textColor=INK, leftIndent=0,
            firstLineIndent=0, spaceBefore=8, spaceAfter=4,
        )]
        toc.dotsMinLevel = 0
        story.append(toc)
    for index, section in enumerate(metadata["sections"], 1):
        story.append(PageBreak())
        story.append(Paragraph(f"SECTION {index:02d}", styles["eyebrow"]))
        heading = Paragraph(rich_text(section["title"]), styles["section"])
        heading.section_key = f"section-{index}"
        story.append(heading)
        if "blocks" in section:
            for block in section["blocks"]:
                story.extend(render_block(block, styles, content_path.parent))
            if section.get("source_note"):
                story.append(Paragraph("<i>" + rich_text(section["source_note"]) + "</i>", styles["small"]))
            continue
        if section.get("callout"):
            story.append(render_callout(section["callout"], styles))
        for paragraph in section.get("paragraphs", []):
            story.append(Paragraph(rich_text(paragraph), styles["body"]))
        if section.get("table"):
            # Dense narrative uses compact vertical padding so a short table
            # need not leave its final row alone on a continuation page.
            padding = 4 if len(section.get("paragraphs", [])) >= 6 else 8
            story.append(render_table(section["table"], styles, padding))
        if section.get("figure"):
            story.append(render_figure(section["figure"], styles, content_path.parent))
        for bullet in section.get("bullets", []):
            story.append(Paragraph("• " + rich_text(bullet), styles["bullet"]))
        if section.get("source_note"):
            story.append(Paragraph("<i>" + rich_text(section["source_note"]) + "</i>", styles["small"]))
    doc.multiBuild(story, maxPasses=5)
    return metadata


def markdown_text(value):
    text = str(value)
    text = re.sub(r'<link\s+[^>]*?href=["\']([^"\']+)["\'][^>]*>(.*?)</link>',
                  lambda match: f"[{match.group(2)}]({match.group(1)})", text,
                  flags=re.S | re.I)
    text = re.sub(r'<font\s+[^>]*?name=["\']DejaVuSansMono[^"\']*["\'][^>]*>(.*?)</font>',
                  lambda match: "`" + match.group(1) + "`", text, flags=re.S | re.I)
    text = re.sub(r"</?font\b[^>]*>", "", text, flags=re.I)
    text = text.replace("<super>", "<sup>").replace("</super>", "</sup>")
    return text


def write_markdown(metadata, output_path):
    """Export an editable source version without changing the supplied facts."""
    lines = ["# " + markdown_text(metadata["title"]), ""]
    if metadata.get("subtitle"):
        lines.extend([markdown_text(metadata["subtitle"]), ""])
    if metadata.get("date"):
        lines.extend([str(metadata["date"]), ""])
    if metadata.get("cover_note"):
        lines.extend(["> " + markdown_text(metadata["cover_note"]), ""])

    def add_table(spec):
        def cell(value):
            return markdown_text(value).replace("|", "\\|").replace("\n", "<br/>")
        lines.append("| " + " | ".join(cell(v) for v in spec["headers"]) + " |")
        lines.append("| " + " | ".join("---" for _ in spec["headers"]) + " |")
        lines.extend("| " + " | ".join(cell(v) for v in row) + " |" for row in spec["rows"])
        lines.append("")

    def add_figure(spec):
        caption = markdown_text(spec.get("caption", "Figure"))
        lines.extend([f"![{plain_text(caption)}]({spec['path']})", "", caption, ""])

    def add_bullets(items):
        lines.extend("- " + markdown_text(v) for v in items)
        lines.append("")

    for index, section in enumerate(metadata["sections"], 1):
        lines.extend([f"## {index}. " + markdown_text(section["title"]), ""])
        if "blocks" in section:
            for block in section["blocks"]:
                kind = block.get("type", "paragraph")
                if kind == "table":
                    add_table(block)
                elif kind == "figure":
                    add_figure(block)
                elif kind == "bullets":
                    add_bullets(block["items"])
                elif kind == "subheading":
                    lines.extend(["### " + markdown_text(block["text"]), ""])
                elif kind == "callout":
                    lines.extend(["> " + markdown_text(block["text"]), ""])
                elif kind == "page_break":
                    lines.extend(["<!-- page break -->", ""])
                else:
                    lines.extend([markdown_text(block["text"]), ""])
        else:
            if section.get("callout"):
                lines.extend(["> " + markdown_text(section["callout"]), ""])
            for paragraph in section.get("paragraphs", []):
                lines.extend([markdown_text(paragraph), ""])
            if section.get("table"):
                add_table(section["table"])
            if section.get("figure"):
                add_figure(section["figure"])
            if section.get("bullets"):
                add_bullets(section["bullets"])
        if section.get("source_note"):
            lines.extend(["*" + markdown_text(section["source_note"]) + "*", ""])
    Path(output_path).write_text("\n".join(lines), encoding="utf-8")


def quality_check(pdf_path, output_dir):
    """Rasterize and audit PDF page boxes and the supplied font coverage."""
    import fitz
    from PIL import Image as PILImage, ImageDraw, ImageFont

    pdf_path, output_dir = Path(pdf_path), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Rebuilding must not leave page-numbered previews from an older layout.
    for stale in output_dir.glob("page_[0-9][0-9].png"):
        stale.unlink()
    document = fitz.open(pdf_path)
    outside_page = []
    replacement_chars = []
    pages = []
    link_count = 0
    for index, page in enumerate(document):
        page_rect = page.rect
        link_count += len(page.get_links())
        for block in page.get_text("dict")["blocks"]:
            if block.get("type") != 0:
                continue
            for line in block["lines"]:
                for span in line["spans"]:
                    x0, y0, x1, y1 = span["bbox"]
                    if x0 < -0.5 or y0 < -0.5 or x1 > page_rect.width + 0.5 or y1 > page_rect.height + 0.5:
                        outside_page.append({"page": index + 1, "text": span["text"], "bbox": span["bbox"]})
                    if "\ufffd" in span["text"] or "\x00" in span["text"]:
                        replacement_chars.append({"page": index + 1, "text": span["text"]})
        pixmap = page.get_pixmap(matrix=fitz.Matrix(0.55, 0.55), alpha=False)
        img = PILImage.frombytes("RGB", [pixmap.width, pixmap.height], pixmap.samples)
        pages.append(img)
        # Per-page previews let a reviewer inspect any dense table or figure.
        preview = page.get_pixmap(matrix=fitz.Matrix(1.5, 1.5), alpha=False)
        preview.save(str(output_dir / f"page_{index + 1:02d}.png"))
    columns = 4
    gap, label_h = 18, 24
    thumb_w, thumb_h = pages[0].size
    rows = math.ceil(len(pages) / columns)
    sheet = PILImage.new("RGB", (gap + columns * (thumb_w + gap),
                                 gap + rows * (thumb_h + label_h + gap)), "#E8EDF1")
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
    for i, img in enumerate(pages):
        x = gap + (i % columns) * (thumb_w + gap)
        y = gap + (i // columns) * (thumb_h + label_h + gap)
        sheet.paste(img, (x, y))
        draw.text((x, y + thumb_h + 5), f"Page {i + 1:02d}", fill="#25354A", font=font)
    contact_path = output_dir / "contact_sheet.jpg"
    sheet.save(contact_path, quality=90)
    summary = {
        "pdf": str(pdf_path.resolve()), "pages": len(document),
        "clickable_links": link_count, "page_bounds_violations": outside_page,
        "replacement_glyphs": replacement_chars,
        "missing_font_glyphs": [{"character": c, "codepoint": f"U+{ord(c):04X}"}
                                for c in sorted(MISSING_GLYPHS)],
        "contact_sheet": str(contact_path.resolve()),
    }
    (output_dir / "layout_qa.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    document.close()
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--content", type=Path, default=ROOT / "report_content.json")
    parser.add_argument("--output", type=Path, default=ROOT / "unified_curve_basis_report.pdf")
    parser.add_argument("--qa-dir", type=Path, default=ROOT / "preview")
    parser.add_argument("--markdown", type=Path, default=ROOT / "report_source.md")
    parser.add_argument("--no-qa", action="store_true", help="Skip page raster previews and bounds/glyph audit")
    args = parser.parse_args()
    metadata = build_pdf(args.content, args.output)
    write_markdown(metadata, args.markdown)
    if args.no_qa:
        print(f"Created {args.output}")
    else:
        qa = quality_check(args.output, args.qa_dir)
        print(json.dumps({k: v for k, v in qa.items() if k != "missing_font_glyphs"}, indent=2))
        if qa["page_bounds_violations"] or qa["replacement_glyphs"] or qa["missing_font_glyphs"]:
            raise SystemExit("Layout QA found bounds or glyph issues; see layout_qa.json")


if __name__ == "__main__":
    main()
