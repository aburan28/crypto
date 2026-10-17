"""Document export only. Experiment arithmetic and statistics are computed natively."""
from pathlib import Path
import sys
import json
import fitz
from playwright.sync_api import sync_playwright

root = Path(sys.argv[1])
charts = json.loads((root / 'charts.json').read_text())
errors = []
pages = []
with sync_playwright() as api:
    browser = api.chromium.launch(executable_path='/usr/bin/chromium', headless=True,
        args=['--no-sandbox','--disable-dev-shm-usage','--disable-crash-reporter'])
    page = browser.new_page(viewport={'width':1600,'height':1100})
    page.route('**/*', lambda route: route.abort())
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.set_content((root / 'report.html').read_text())
    assert page.locator('svg .node').count() == len(charts)*46
    assert page.locator('svg .edge').count() == len(charts)*90
    assert not (root / 'RUN_REPORT.pdf').exists()
    page.pdf(path=str(root / 'RUN_REPORT.pdf'),print_background=True,prefer_css_page_size=True)
    # A short results companion selects the actual observations and two key
    # diagnostic figures from exactly the same report source; no recomputation.
    page.evaluate("""() => {
      document.querySelector('h1').textContent='Photographed-curve restart: results and limits';
      document.querySelectorAll('section.page').forEach((section,index)=>{
        if (![4,5,7].includes(index)) section.remove();
      });
    }""")
    assert not (root / 'RESULTS_REPORT.pdf').exists()
    page.pdf(path=str(root / 'RESULTS_REPORT.pdf'),print_background=True,prefer_css_page_size=True)
    for chart in charts:
        directory = root / chart['name']
        image_page = browser.new_page(viewport={'width':2600,'height':1400},device_scale_factor=150/96)
        image_page.route('**/*', lambda route: route.abort())
        svg = (directory / 'volcano.svg').read_text()
        image_page.set_content('<style>body{margin:0}svg{display:block;width:2440px;height:auto}</style>'+svg)
        image_page.locator('svg').screenshot(path=str(directory / 'volcano.png'))
        image_page.close()
    browser.close()
for filename in ['RUN_REPORT.pdf','RESULTS_REPORT.pdf']:
    document = fitz.open(root / filename)
    for i, page in enumerate(document):
        name = filename.removesuffix('.pdf')+f'_page_{i+1}.png'
        page.get_pixmap(matrix=fitz.Matrix(1.2,1.2)).save(root / name)
        pages.append({'document':filename,'page':i+1,'text_chars':len(page.get_text()),'image':name})
assert not errors
(root / 'render_audit.json').write_text(json.dumps({'passed':True,'pages':pages,
    'figures':len(charts),'nodes_per_figure':46,'directed_maps_per_figure':90,
    'png_scale':150/96,'javascript_errors':errors,'human_page_inspection_required':True},indent=2)+'\n')
print(json.dumps({'pages':pages,'figures':len(charts)}))
