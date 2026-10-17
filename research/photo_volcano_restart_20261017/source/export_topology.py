"""Document rendering/inspection only; no cryptographic computation."""
from pathlib import Path
import sys,json
import fitz
from playwright.sync_api import sync_playwright
report=Path(sys.argv[1]);errors=[]
with sync_playwright() as api:
 browser=api.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--disable-dev-shm-usage','--disable-crash-reporter'])
 page=browser.new_page(viewport={'width':1600,'height':1100},device_scale_factor=150/96);page.route('**/*',lambda r:r.abort());page.on('pageerror',lambda e:errors.append(str(e)));page.set_content((report/'report.html').read_text());assert page.locator('svg .node').count()==46;assert page.locator('svg .edge').count()==90;page.pdf(path=str(report/'REPORT.pdf'),print_background=True,prefer_css_page_size=True);page.locator('svg').screenshot(path=str(report/'volcano.png'));browser.close()
doc=fitz.open(report/'REPORT.pdf');pages=[]
for i,page in enumerate(doc):
 path=report/f'page_{i+1}.png';page.get_pixmap(matrix=fitz.Matrix(1.25,1.25)).save(path);pages.append({'page':i+1,'text_chars':len(page.get_text()),'image':path.name})
assert not errors
(report/'render_audit.json').write_text(json.dumps({'passed':True,'nodes':46,'directed_maps':90,'pages':pages,'png_scale':150/96,'javascript_errors':errors},indent=2)+'\n');print(json.dumps(pages))
