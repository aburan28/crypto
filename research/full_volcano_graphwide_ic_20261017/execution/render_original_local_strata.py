"""Present native local-conductor certificates; no mathematical aggregation."""
import json
import sys
from pathlib import Path
import typst

data,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
certs=[(label,json.loads((data/f'{label}.json').read_text())) for label in 'ABCDE']
(out/'plot-data.json').write_text(json.dumps({label:{k:v for k,v in c.items() if k!='rows'} for label,c in certs},indent=2)+'\n')
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="540" viewBox="0 0 1100 540"><rect width="1100" height="540" fill="white"/><g font-family="sans-serif" fill="#18364c">',
     '<text x="20" y="30" font-size="23" font-weight="bold">Original A–E: locally certified prefix conductors</text>',
     '<text x="20" y="55" font-size="14">Blue: exact local conductor. Grey: unknown. Complete class populations are not established by these prefixes.</text>']
for i,(label,c) in enumerate(certs):
    y=85+85*i
    slug=c['plan']['source_identity']['slug']
    svg += [f'<text x="20" y="{y}" font-size="15" font-weight="bold">{label}: {slug}</text>',
            f'<rect x="20" y="{y+12}" width="{c["vertices_exact"]*.34}" height="27" fill="#256f9a"/>',
            f'<rect x="{20+c["vertices_exact"]*.34}" y="{y+12}" width="{c["vertices_unknown"]*.34}" height="27" fill="#cdd8de"/>',
            f'<text x="690" y="{y+31}" font-size="14">{c["vertices_exact"]} exact; {c["vertices_unknown"]} unknown; {c["strata_exact_observed"]} observed strata</text>']
svg+=['<text x="20" y="526" font-size="13">Same ordinary-volcano theorem and map arithmetic; calibrated on complete graphs. No independent arithmetic or formal proof claim.</text></g></svg>']
(out/'coverage.svg').write_text('\n'.join(svg)+'\n')
rows=''
for label,c in certs:
    cells=[label,c['plan']['input']['p'],c['rows'][0]['conductor'],c['vertices'],c['vertices_exact'],c['strata_exact_observed']]
    rows+=', '.join(f'[{v}]' for v in cells)+',\n'
source=r'''#set page(paper:"a4",margin:13mm)
#set text(size:9pt)
= Original A--E: local End(E) conductor certificates
17 October 2026. Input certification only; no new ECDLP timing run.
#image("coverage.svg",width:100%)
#table(columns:(.3fr,1.8fr,1.2fr,.6fr,.6fr,.6fr),inset:3pt,
 [Case], [Prime field p], [Source conductor], [Prefix vertices], [Exact vertices], [Observed strata],
''' + rows + r''')

The conductor is the requested index of End(E) in the maximal CM order.
These source values are locally determined from complete rational-kernel
counts and the standard ordinary-volcano floor-distance correspondence;
they are not merely copied Frobenius conductors. For each prime, distance
to known floors bounds height above; distance to known floors or unexpanded
frontier bounds it below. Only equal finite bounds produce exact labels.
Unknown vertices are never declared floors or omitted from coverage.

The same calculation passed all available complete-truth bound/label checks
on the three complete additional graphs at visibility cutoffs 0, 1, 17, 1024
and all. The generic certify mode reproduces the earlier C certificate byte
for byte. It shares isogeny arithmetic and the standard theorem; this is not
an independent arithmetic implementation or a new mathematical theorem.

All five frozen prefix DBs are read immutably. Every output binds its DB
SHA-256 and existing map-replay receipt. The original databases and global
stratum status are unchanged. Actual local certificate rows and per-conductor
counts are in the accompanying archive. No prefix is labeled a complete class.
Original full A--E total 7,833,658,967 expected vertices and 515 strata;
historical B has a separate complete delivery. A/C/D/E remain incomplete.

The largest prime subgroup and embedding degree are separate properties.
In particular C has largest prime order 7 despite its 61-bit field and
source conductor 967,458,816. These local labels do not establish a speed
advantage, a vulnerability caused by conductor, a full benchmark matrix or
finite-field solver support for all models. No deployed-security claim.
'''
(out/'report.typ').write_text(source)
typst.compile(str(out/'report.typ'),output=str(out/'report.pdf'))
