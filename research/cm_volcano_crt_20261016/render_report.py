"""Render recorded Markdown/SVG evidence only; no algorithms or measurements."""
from pathlib import Path
import re
from xml.sax.saxutils import escape
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from svglib.svglib import svg2rlg
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

pdfmetrics.registerFont(TTFont('ReportSans', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('ReportBold', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('ReportMono', '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))

root = Path(__file__).resolve().parent
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='BodyReport', fontName='ReportSans', fontSize=9, leading=13, spaceAfter=7))
styles.add(ParagraphStyle(name='CellReport', fontName='ReportSans', fontSize=7.5, leading=10.5, spaceAfter=0))
styles['Title'].fontName = 'ReportBold'
styles['Title'].fontSize = 18
styles['Heading2'].fontName = 'ReportBold'
styles['Title'].textColor = colors.HexColor('#17334d')
styles['Heading2'].textColor = colors.HexColor('#17334d')
styles['Heading2'].spaceBefore = 13
styles['Heading2'].spaceAfter = 7

def inline(s):
    s = escape(s)
    s = re.sub(r'`([^`]+)`', r'<font name="ReportMono">\1</font>', s)
    s = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<font color="#24587a">\1</font>', s)
    return s

def paragraph(s, cell=False):
    return Paragraph(inline(s), styles['CellReport' if cell else 'BodyReport'])

lines = (root / 'REPORT.md').read_text().splitlines()
story = []
i = 0
while i < len(lines):
    line = lines[i].strip()
    if not line:
        i += 1
        continue
    if line.startswith('# '):
        story.append(Paragraph(inline(line[2:]), styles['Title']))
    elif line.startswith('## '):
        story.append(Paragraph(inline(line[3:]), styles['Heading2']))
    elif line.startswith('!['):
        asset = re.search(r'\]\(([^)]+)\)', line).group(1)
        drawing = svg2rlg(str(root / asset))
        scale = 510 / drawing.width
        drawing.scale(scale, scale)
        drawing.width *= scale
        drawing.height *= scale
        story.extend([Spacer(1,5), drawing, Spacer(1,7)])
    elif line.startswith('|'):
        rows = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            cells = lines[i].strip().strip('|').split('|')
            if not all(re.fullmatch(r'\s*:?-+:?\s*', x) for x in cells):
                rows.append([paragraph(x.strip(), True) for x in cells])
            i += 1
        n = len(rows[0])
        widths = [95,160,255] if n == 3 else [510/n]*n
        table = Table(rows, colWidths=widths, repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),colors.HexColor('#dceaf4')),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,colors.HexColor('#f4f7fa')]),
            ('VALIGN',(0,0),(-1,-1),'TOP'),
            ('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),
            ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
            ('LINEBELOW',(0,0),(-1,0),0.7,colors.HexColor('#416c8d')),
        ]))
        story.extend([table,Spacer(1,8)])
        continue
    else:
        parts = [line]
        while i+1 < len(lines) and lines[i+1].strip() and not lines[i+1].lstrip().startswith(('#','|','![','- ')):
            i += 1
            parts.append(lines[i].strip())
        story.append(paragraph(' '.join(parts)))
    i += 1

def footer(canvas, doc):
    canvas.setFont('ReportSans',7)
    canvas.setFillColor(colors.HexColor('#536878'))
    canvas.drawString(42,25,'CM volcano CRT - 2026-10-16 - exact correctness evidence; no performance claim')
    canvas.drawRightString(553,25,str(doc.page))

SimpleDocTemplate(str(root / 'REPORT.pdf'),pagesize=(595,842),leftMargin=42,rightMargin=42,
    topMargin=38,bottomMargin=43,title='CM volcano modular polynomials and explicit CRT',author='crypto repository').build(
    story,onFirstPage=footer,onLaterPages=footer)
