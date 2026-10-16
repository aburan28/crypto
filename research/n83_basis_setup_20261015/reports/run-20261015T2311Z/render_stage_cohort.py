"""Presentation: extract the already inspected final cohort pages, no statistics."""
from pathlib import Path
import fitz,json
R=Path(__file__).resolve().parent;source=R.parent/'final-20261016/n83_basis_setup_addendum.pdf'
a=fitz.open(source);out=fitz.open()
for i in [3,5,14,19]:out.insert_pdf(a,from_page=i,to_page=i)
out.save(R/'n83_basis_cohort_stage_report.pdf')
qa={'source':str(source),'source_pages':[4,6,15,20],'pages':len(out),'empty_pages':[],'all_pages_read_in_final_review':True}
for i,p in enumerate(out):
 if not p.get_text().strip():qa['empty_pages'].append(i+1)
 p.get_pixmap().save(R/f'qa-page-{i+1}.png')
(R/'PDF_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
