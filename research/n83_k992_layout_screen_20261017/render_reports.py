#!/usr/bin/env python3
"""Render deterministic run and result PDFs from the sealed replay."""

import json
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / 'replay.json').read_text())
assert DATA['status'] == 'PASS_INDEPENDENT_ARITHMETIC_REPLAY'
WIDTH, HEIGHT = letter
INK = colors.HexColor('#173042')
MUTED = colors.HexColor('#536875')
ACCENT = colors.HexColor('#217a9a')
LINE = colors.HexColor('#cbd8df')


def begin(filename, title, subtitle):
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


def lines(pdf, y, text_lines, leading=14):
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica', 9)
    for line in text_lines:
        pdf.drawString(36, y, line)
        y -= leading
    return y


def footer(pdf, label):
    pdf.setStrokeColor(LINE)
    pdf.line(36, 44, 576, 44)
    pdf.setFillColor(MUTED)
    pdf.setFont('Helvetica', 8)
    pdf.drawString(36, 30, label)
    pdf.drawRightString(576, 30, '2026-10-17  |  1 / 1')
    pdf.showPage()
    pdf.save()


def run_report():
    pdf = begin('RUN.pdf', 'N83 K992 layout control - run record',
                'Four fresh processes on one previously exposed public point')
    y = heading(pdf, 663, 'Instance and protocol')
    y = lines(pdf, y, [
        'Curve: GF(2^83), y^2 + xy = x^3 + x^2 + 1; prime subgroup order 8569786107849059.',
        'Fixed point: research/n83_cold_factor_base_20261016/targets.jsonl.',
        'Arms: K928/root11, K960/root11, K992/root11, K992/root22.',
        'Each process rebuilt the base, dynamic 75%-load root index, and factor-base logs.',
        'Symmetric S3, sparse canonical map, four rank threads, one target thread.',
        'CPU: shared Apple M4 Pro, L0; timing is exploratory, with no isolation gate.',
    ])
    y -= 13
    y = heading(pdf, y, 'Source and immutable inputs')
    y = lines(pdf, y, [
        'IC source SHA-256:',
        DATA['source_sha256'],
        'Executed binary SHA-256:',
        DATA['binary_sha256'],
        'Exposed point SHA-256:',
        DATA['target_sha256'],
    ], leading=15)
    y -= 12
    y = heading(pdf, y, 'Verification and custody')
    y = lines(pdf, y, [
        'Independent replay rebuilt every base point and summed each four-point relation.',
        'It checked the recovered scalar and the exclusive target-online phase sum.',
        'All four rows: PASS_INDEPENDENT_ARITHMETIC_REPLAY.',
        'The 28-member raw archive was checked against the versioned manifest.',
        'Raw archive SHA-256:',
        '02993c7b3b7b7bf85af28bf999acfc48b41df2d61afa9408dfde3d9f440d5f55',
    ], leading=15)
    footer(pdf, 'Research record: research/n83_k992_layout_screen_20261017')


def result_report():
    pdf = begin('RESULT.pdf', 'N83 K992 layout control - results',
                'Complete cold process and one-target online lookup, same exposed point')
    y = heading(pdf, 663, 'Measured intervals')
    headers = ['Base', 'Usable pts', 'Index s', 'Rank s', 'Cold s', 'Online ms']
    x = [42, 143, 251, 322, 393, 470]
    pdf.setFillColor(ACCENT)
    pdf.rect(36, y - 11, 540, 25, fill=1, stroke=0)
    pdf.setFillColor(colors.white)
    pdf.setFont('Helvetica-Bold', 8)
    for start, label in zip(x, headers):
        pdf.drawString(start, y, label)
    y -= 27
    pdf.setFont('Helvetica', 8)
    for index, row in enumerate(DATA['rows']):
        if index % 2 == 0:
            pdf.setFillColor(colors.HexColor('#edf4f7'))
            pdf.rect(36, y - 6, 540, 23, fill=1, stroke=0)
        pdf.setFillColor(INK)
        values = [f"K{row['k']}/" + ('root22' if row['layout'] == 'root22_dynamic' else 'root11'),
                  f"{row['actual_usable_points']:,}",
                  f"{row['index_build_ms']/1000:.3f}",
                  f"{row['rank_stage_ms']/1000:.3f}",
                  f"{row['cold_process_ms']/1000:.3f}",
                  f"{row['target_online_ms']:.3f}"]
        for start, value in zip(x, values):
            pdf.drawString(start, y + 2, value)
        y -= 24
    y -= 10
    pdf.setFillColor(INK)
    pdf.setFont('Helvetica', 9)
    pdf.drawString(36, y, 'K992/root11: 75.841 s cold; root22: 76.330 s. The layouts are tied here.')
    y -= 18
    pdf.drawString(36, y, 'K960 uses 680,614 target probes; K992 uses 117,653 on this point.')
    image_y = 148
    pdf.drawImage(str(HERE / 'screen.png'), 36, image_y, width=540,
                  height=228, preserveAspectRatio=True, anchor='c')
    pdf.setFillColor(MUTED)
    pdf.setFont('Helvetica', 8)
    pdf.drawString(36, 119, 'L0 local wall times vary: an earlier K992/root11 run took 95.223 s.')
    pdf.drawString(36, 106, 'The sealed ten-point Tier 1 K992/K1024/K1107 panel is the fresh-point gate.')
    pdf.drawString(36, 93, 'Every row passed independent base, relation, scalar and phase replay.')
    footer(pdf, 'Research record: research/n83_k992_layout_screen_20261017')


run_report()
result_report()
print(json.dumps({'status': 'PASS_RENDERED', 'pdfs': ['RUN.pdf', 'RESULT.pdf']},
                 sort_keys=True))
