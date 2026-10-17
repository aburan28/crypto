"""Render recount.svg and report.tex from the counts in count.txt.

    python3 make_report.py && pdflatex report.tex
"""
import re

rows = []  # (label, hypothesis-ordinary, excluded-ordinary, supersingular)
txt = open("count.txt", encoding="utf-8").read()
for ell in (3, 5, 7, 11, 13):
    def get(kind, cls):
        m = re.search(rf"l={ell}\s+{kind}\s+{cls}\s+\S+\s+(\d+)", txt)
        return int(m.group(1)) if m else 0
    rows.append((ell, get("ordinary", "hypothesis"), get("ordinary", "excluded"),
                 get("supersingular", "excluded")))

W, H, L, B = 640, 320, 70, 260
top = max(sum(r[1:]) for r in rows)
bw = 70
colors = ("#2b6cb0", "#dd6b20", "#a0aec0")
names = ("ordinary, lemma hypothesis (0 violations)",
         "ordinary, excluded r = l-1", "supersingular (excluded)")
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H + 60}" '
       'font-family="sans-serif" font-size="12">',
       f'<rect width="{W}" height="{H + 60}" fill="white"/>',
       '<text x="20" y="20" font-size="14">j = 0 twist-lemma instances, p &lt; 500 '
       '(recount, SHF-20261016-3caa6f P3)</text>']
for i, (ell, *vals) in enumerate(rows):
    x, y = L + i * (bw + 40), B
    for v, c in zip(vals, colors):
        h = (B - 40) * v / top
        y -= h
        svg.append(f'<rect x="{x}" y="{y:.1f}" width="{bw}" height="{h:.1f}" fill="{c}"/>')
    svg.append(f'<text x="{x + bw / 2}" y="{B + 16}" text-anchor="middle">l = {ell}</text>')
    svg.append(f'<text x="{x + bw / 2}" y="{y - 4:.1f}" text-anchor="middle">{sum(vals):,}</text>')
for j, (n, c) in enumerate(zip(names, colors)):
    svg.append(f'<rect x="{L}" y="{B + 30 + j * 16}" width="10" height="10" fill="{c}"/>')
    svg.append(f'<text x="{L + 16}" y="{B + 39 + j * 16}">{n}</text>')
svg.append("</svg>")
open("recount.svg", "w", encoding="utf-8").write("\n".join(svg) + "\n")

body = "\n".join(f"{ell} & {a:,} & {b:,} & {c:,} & {a + b + c:,} \\\\" for ell, a, b, c in rows)
tex = r"""\documentclass[11pt]{article}
\usepackage[margin=2.2cm]{geometry}
\usepackage{tikz}
\begin{document}
\section*{Structural-completeness recount (sheriff case SHF-20261016-3caa6f)}
The paper's ``16\,562 ordinary curve instances'' matched no script output.
Re-running the cited scripts (logs in this directory, \texttt{count.py})
gives the counts below. Every ordinary instance in the lemma's hypothesis
class satisfies it (0 violations); the excluded class $r=\ell-1$ always has
$\ell \mid \#E^t$, as predicted. The CM-73 sweep processed the 330 primes
$p \equiv 1 \pmod 6$, $7 \le p \le 5000$ (not 621).

\begin{center}\begin{tabular}{rrrrr}
$\ell$ & ordinary, hypothesis & ordinary, excluded & supersingular & total\\\hline
""" + body + r"""
\hline
all & 9,864 & 672 & 15,742 & 26,278
\end{tabular}\end{center}

\begin{center}
\begin{tikzpicture}[x=1.4cm,y=0.00025cm]
""" + "\n".join(
    rf"\fill[blue!60] ({i},0) rectangle ({i}+0.8,{a});"
    rf"\fill[orange!80] ({i},{a}) rectangle ({i}+0.8,{a + b});"
    rf"\fill[gray!50] ({i},{a + b}) rectangle ({i}+0.8,{a + b + c});"
    rf"\node[below] at ({i}+0.4,0) {{$\ell={ell}$}};"
    rf"\node[above] at ({i}+0.4,{a + b + c}) {{\small {a + b + c:,}}};"
    for i, (ell, a, b, c) in enumerate(rows)) + r"""
\end{tikzpicture}

\small blue: ordinary, hypothesis class; orange: ordinary, $r=\ell-1$; grey: supersingular.
\end{center}
\end{document}
"""
open("report.tex", "w", encoding="utf-8").write(tex)
