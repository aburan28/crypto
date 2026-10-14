"""Independent PDF layout/preservation inspection; no research computation."""
from pathlib import Path
import hashlib
import json
import fitz
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parent
ADDENDUM = Path('/workspace/report_output/n83_faster_rho_evaluation.pdf')
COMBINED = Path('/workspace/report_output/unified_curve_basis_report_faster_rho.pdf')
OLD = Path('/workspace/report_output/unified_curve_basis_report_optimized_rho.pdf')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def inspect(path):
    doc = fitz.open(path)
    bounds = []
    glyphs = []
    links = []
    texts = []
    for i, page in enumerate(doc):
        texts.append(page.get_text())
        for link in page.get_links():
            links.append({'page': i + 1, 'uri': link.get('uri'), 'kind': link['kind']})
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                for span in line['spans']:
                    box = fitz.Rect(span['bbox'])
                    if box.x0 < -0.5 or box.y0 < -0.5 or box.x1 > page.rect.width + 0.5 or box.y1 > page.rect.height + 0.5:
                        bounds.append({'page': i + 1, 'bbox': list(box), 'text': span['text']})
                    if '\ufffd' in span['text'] or '\x00' in span['text']:
                        glyphs.append({'page': i + 1, 'text': span['text']})
    return doc, {'file': str(path), 'sha256': sha(path.read_bytes()), 'pages': len(doc),
                 'clickable_links': len(links), 'bounds_violations': bounds,
                 'replacement_or_null_glyphs': glyphs, 'links': links}, texts

add, a, addtexts = inspect(ADDENDUM)
combined, c, _ = inspect(COMBINED)
old, o, _ = inspect(OLD)
preserved = []
for i in range(len(old)):
    historical = old[i]
    appended = combined[len(add) + i]
    hp = historical.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
    cp = appended.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
    preserved.append({'old_page': i + 1, 'combined_page': len(add) + i + 1,
                      'text_equal': historical.get_text() == appended.get_text(),
                      'size_equal': (hp.width, hp.height, hp.n) == (cp.width, cp.height, cp.n),
                      'pixels_equal_at_144dpi': hp.samples == cp.samples,
                      'old_pixels_sha256': sha(hp.samples), 'combined_pixels_sha256': sha(cp.samples)})

thumbs = []
for i, page in enumerate(add):
    pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
    pix.save(str(ROOT / f'page_{i+1:02d}.png'))
    im = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
    im.thumbnail((430, 600))
    thumbs.append(im)
sheet = Image.new('RGB', (3 * 450, ((len(thumbs) + 2) // 3) * 630), 'white')
draw = ImageDraw.Draw(sheet)
for i, im in enumerate(thumbs):
    x, y = (i % 3) * 450, (i // 3) * 630
    draw.text((x + 8, y + 4), f'Addendum page {i+1}', fill='black')
    sheet.paste(im, (x + 8, y + 25))
sheet.save(ROOT / 'contact_sheet.jpg', quality=92)

fulltext = '\n\n'.join(f'PAGE {i+1}\n{text}' for i, text in enumerate(addtexts))
(ROOT / 'addendum_text.txt').write_text(fulltext)
expected_strings = ['1.069651', '1.043178', '1.089238', '2.987008', '2.780545',
                    '2.923128', '2.903055', '0.962508', '1.046541', '139',
                    '256', '64', '60', 'manual', 'Lean']
result = {'inspection_kind': 'pdf_layout_and_historical_pixel_preservation_only',
          'research_results_recomputed': False, 'addendum': a, 'combined': c,
          'historical': o, 'page_count_sum_matches': len(combined) == len(add) + len(old),
          'expected_report_strings_present': {s: s in fulltext for s in expected_strings},
          'historical_page_preservation_at_144dpi': preserved,
          'all_60_historical_pages_preserved': len(old) == 60 and all(p['text_equal'] and p['size_equal'] and p['pixels_equal_at_144dpi'] for p in preserved),
          'visual_review_status': 'pending_model_inspection_of_every_new_page'}
(ROOT / 'qa.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'addendum_pages': len(add), 'combined_pages': len(combined), 'historical_pages': len(old),
                  'all_old_pages_pixel_identical': result['all_60_historical_pages_preserved'],
                  'bounds_violations': len(a['bounds_violations']),
                  'glyph_issues': len(a['replacement_or_null_glyphs']),
                  'required_text_present': all(result['expected_report_strings_present'].values())}, indent=2))
