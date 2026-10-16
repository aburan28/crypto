"""Presentation QA: text and pixel preservation with independent renderers."""
import fitz,hashlib,json,subprocess
from pathlib import Path
R=Path(__file__).resolve().parent
old=Path('/workspace/report_output/unified_curve_basis_report_matched_ic_rho.pdf');new=R/'unified_curve_basis_report_basis_setup.pdf';add=R/'n83_basis_setup_addendum.pdf'
assert hashlib.sha256(old.read_bytes()).hexdigest()=='af881641988a98bec3177e286694233f7b03731e82b6517a51ea78fc5e19793a'
a,b,c=fitz.open(old),fitz.open(new),fitz.open(add);assert len(a)==87 and len(b)==87+len(c)
rows=[]
for i in range(len(b)):
 source=a if i<87 else c;sourcepath=old if i<87 else add;j=i if i<87 else i-87
 left,right=source[j],b[i];t=left.get_text()==right.get_text();x=left.get_pixmap(matrix=fitz.Matrix(1,1));y=right.get_pixmap(matrix=fitz.Matrix(1,1));p=x.samples==y.samples and x.width==y.width and x.height==y.height
 def independent(path,page):return subprocess.check_output(['pdftoppm','-f',str(page+1),'-l',str(page+1),'-r','144','-singlefile',str(path)],stderr=subprocess.DEVNULL)
 u,v=independent(sourcepath,j),independent(new,i);z=u==v
 rows.append({'comprehensive_page':i+1,'source_page':j+1,'source':'historical87' if i<87 else 'new_addendum','text_equal':t,'fitz72_pixels_equal':p,'independent_poppler144_pixels_equal':z,'pixel144_sha256':hashlib.sha256(u).hexdigest()})
 assert t and p and z,(i+1,t,p,z)
 if i%20==0:print('Preserved pages through',i+1,flush=True)
receipt={'status':'PASS','historical_pages':87,'new_pages':len(c),'comprehensive_pages':len(b),'old_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'new_sha256':hashlib.sha256(new.read_bytes()).hexdigest(),'all_pages_text_and_pixels72_and_independent144_preserved':True,'pages':rows}
(R/'PDF_PRESERVATION.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS preservation',len(b),'pages',flush=True)
