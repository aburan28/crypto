"""Document rendering only. Does not compute experiment results."""
from pathlib import Path
import json
import sys
import fitz
from playwright.sync_api import sync_playwright

root = Path(sys.argv[1])
for name in ['REPORT.pdf', 'volcano.png', 'render_audit.json']:
    assert not (root / name).exists(), f'refuse overwrite {name}'
manifest = json.loads((root / 'render_manifest.json').read_text())
errors = []
with sync_playwright() as api:
    browser = api.chromium.launch(executable_path='/usr/bin/chromium', headless=True,
        args=['--no-sandbox', '--disable-dev-shm-usage', '--disable-crash-reporter'])
    page = browser.new_page(viewport={'width':1600, 'height':1000}, device_scale_factor=150/96)
    page.route('**/*', lambda route: route.abort())
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.set_content((root / 'index.html').read_text())
    assert page.locator('svg .node').count() == manifest['node_count']
    assert page.locator('svg .edge').count() == manifest['directed_edge_count']
    page.locator('svg').screenshot(path=str(root / 'volcano.png'))
    page.pdf(path=str(root / 'REPORT.pdf'), print_background=True, prefer_css_page_size=True)
    browser.close()
pages = []
for i, page in enumerate(fitz.open(root / 'REPORT.pdf')):
    filename = f'page_{i+1}.png'
    page.get_pixmap(matrix=fitz.Matrix(1.2,1.2)).save(root / filename)
    pages.append({'page':i+1,'text_chars':len(page.get_text()),'image':filename})
assert not errors
(root / 'render_audit.json').write_text(json.dumps({'passed':True,'pages':pages,
    'nodes':manifest['node_count'],'directed_maps':manifest['directed_edge_count'],
    'png_scale':150/96,'javascript_errors':errors,'human_page_inspection_required':True},indent=2)+'\n')
print(json.dumps(pages))
