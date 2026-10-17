#!/usr/bin/env python3
"""Render K984/K1000/K1016 run and results PDFs from replayed data."""

import json
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / 'replay.json').read_text())
assert DATA['status'] == 'PASS_INDEPENDENT_ARITHMETIC_REPLAY'
ROWS = sorted(DATA['rows'], key=lambda row: row['k'])
INK = colors.HexColor('#173042')
MUTED = colors.HexColor('#536875')
ACCENT = colors.HexColor('#217a9a')
LINE = colors.HexColor('#cbd8df')


def page(filename, title, subtitle):
    pdf = canvas.Canvas(str(HERE / filename), pagesize=letter,
                        pageCompression=1, invariant=1)
    pdf.setTitle(title)
    pdf.setAuthor('N83 factor-base experiment')
    pdf.setFillColor(ACCENT)
    pdf.rect(36, 752, 540, 5, fill=1, stroke=0)
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica-Bold', 17)
    pdf.drawString(36, 718, title)
    pdf.setFillColor(MUTED)
    pdf.setFont('Helvetica', 9)
    pdf.drawString(36, 699, subtitle)
    pdf.setStrokeColor(LINE)
    pdf.line(36, 687, 576, 687)
    return pdf


def heading(pdf, y, title):
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica-Bold', 11)
    pdf.drawString(36, y, title)
    return y - 18


def text(pdf, y, content):
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica', 9)
    for line in content:
        pdf.drawString(36, y, line)
        y -= 15
    return y


def finish(pdf):
    pdf.setStrokeColor(LINE)
    pdf.line(36, 44, 576, 44)
    pdf.setFillColor(MUTED)
    pdf.setFont('Helvetica', 8)
    pdf.drawString(36, 30, 'Research record: research/n83_k1000_extension_20261017')
    pdf.drawRightString(576, 30, '2026-10-17  |  1 / 1')
    pdf.showPage()
    pdf.save()


pdf = page('RUN.pdf', 'N83 adjacent-base screen - run record',
           'Three fresh processes on one previously exposed public point')
y = heading(pdf, 663, 'Instance and method')
y = text(pdf, y, [
    'Curve: GF(2^83), y^2 + xy = x^3 + x^2 + 1; prime subgroup order 8569786107849059.',
    'Arms: K984/root11, K1000/root11 and K1016/root11.',
    'Each process rebuilt its base, dynamic 75%-load root index and factor-base logs.',
    'Symmetric S3, sparse canonical map, four rank threads, one target thread.',
    'CPU: shared Apple M4 Pro, L0 exploratory; one previously exposed public point.',
])
y -= 16
y = heading(pdf, y, 'Frozen source and input hashes')
y = text(pdf, y, ['IC source SHA-256:', DATA['source_sha256'],
                  'Executed binary SHA-256:', DATA['binary_sha256'],
                  'Exposed point SHA-256:', DATA['target_sha256']])
y -= 16
y = heading(pdf, y, 'Replay and archive')
y = text(pdf, y, [
    'All three full solves reached rank K and recovered the fixed scalar.',
    'Independent replay rebuilt every base and summed every four-point relation.',
    'It checked scalar replay and the exclusive target-online phase sums.',
    'All rows: PASS_INDEPENDENT_ARITHMETIC_REPLAY.',
    'The 23-member raw archive matched its versioned file manifest.',
    'Raw archive SHA-256:',
    '1ebe3cc221b03ce2003078806f4dc2d8fc6f8c26a7976b4368386f62daced5c3',
])
finish(pdf)

pdf = page('RESULT.pdf', 'N83 adjacent-base screen - results',
           'K984, K1000 and K1016 cold process and target-online intervals')
y = heading(pdf, 663, 'Measured intervals')
headers = ('Base', 'Usable pts', 'Index s', 'Rank s', 'Cold s', 'Online ms')
x = (42, 143, 251, 322, 393, 470)
pdf.setFillColor(ACCENT)
pdf.rect(36, y - 11, 540, 25, fill=1, stroke=0)
pdf.setFillColor(colors.white)
pdf.setFont('Helvetica-Bold', 8)
for left, label in zip(x, headers):
    pdf.drawString(left, y, label)
y -= 27
pdf.setFont('Helvetica', 8)
for index, row in enumerate(ROWS):
    if index % 2 == 0:
        pdf.setFillColor(colors.HexColor('#edf4f7'))
        pdf.rect(36, y - 6, 540, 23, fill=1, stroke=0)
    pdf.setFillColor(INK)
    values = (f"K{row['k']}/root11", f"{row['actual_usable_points']:,}",
              f"{row['index_build_ms']/1000:.3f}",
              f"{row['rank_stage_ms']/1000:.3f}",
              f"{row['cold_process_ms']/1000:.3f}",
              f"{row['target_online_ms']:.3f}")
    for left, value in zip(x, values):
        pdf.drawString(left, y + 2, value)
    y -= 24
pdf.setFillColor(INK)
pdf.setFont('Helvetica', 9)
pdf.drawString(36, y - 12, 'K1016 is shortest here at 94.452 s; K992 earlier was 75.841 s.')
pdf.drawString(36, y - 30, 'The K992 reference is from a different local interval, not a paired estimate.')
pdf.drawImage(str(HERE / 'screen.png'), 36, 148, width=540,
              height=228, preserveAspectRatio=True, anchor='c')
pdf.setFillColor(MUTED)
pdf.setFont('Helvetica', 8)
pdf.drawString(36, 119, 'All three controls passed independent base, relation, scalar and phase replay.')
pdf.drawString(36, 106, 'L0 exposed-point data is diagnostic; the ten-point Tier 1 panel is the selection gate.')
finish(pdf)

print(json.dumps({'status': 'PASS_RENDERED', 'pdfs': ['RUN.pdf', 'RESULT.pdf']},
                 sort_keys=True))
