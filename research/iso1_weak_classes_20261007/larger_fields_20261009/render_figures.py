"""Render frozen GP measurements and documentary curve IDs; no curve computation."""
import json, textwrap
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/'curve_records.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
colors={3:'#235b9d',5:'#16826c',7:'#9a5b17'}
fig,(ax,work)=plt.subplots(2,1,figsize=(7.1,6.2),gridspec_kw={'height_ratios':[1,1.15]},layout='constrained')
for n in (3,5,7):
    rows=[r for r in DATA['records'] if r['odd_degree']==n]
    ax.scatter([float(r['log2_field_cardinality']) for r in rows],[n]*len(rows),s=65,c=colors[n],label=f'n = {n}: {len(rows)} paired counts',zorder=3)
ax.set(xlim=(20,267),yticks=[3,5,7],ylim=(2.2,7.8),xlabel='log2(field cardinality), bits',ylabel='Odd extension degree n',title='Verified norm-one curves over F_(p^(2n))')
ax.grid(axis='x',alpha=.23);ax.legend(loc='center',bbox_to_anchor=(.5,.68),ncol=3,fontsize=8,columnspacing=.8,frameon=False)
lookup={}
for line in (ROOT/'field_sizes.txt').read_text().splitlines():
    p,n,q,bits,c,logc=line.split('|');lookup[(int(p),int(n))]=float(logc)
for n in (3,5,7):
    xs=[13,257,1009]
    work.plot(range(3),[lookup[(p,n)] for p in xs],marker='o',lw=2,color=colors[n],label=f'n = {n}')
work.set(xticks=range(3),xticklabels=['p = 13','p = 257','p = 1009'],ylabel='log10(point-count calls)',ylim=(0,37),title='Exact orbit work for a complete parameter census')
work.grid(alpha=.23);work.legend(fontsize=8,loc='upper left')
fig.suptitle('Larger fields verified; complete-census work remains explicit',fontsize=13)
fig.savefig(ROOT/'field_expansion.svg');fig.savefig(ROOT/'field_expansion.pdf');plt.close(fig)
pair=next(r for r in DATA['records'] if r['characteristic']==13 and r['odd_degree']==5)
fig,ax=plt.subplots(figsize=(7.1,3.1));fig.subplots_adjust(left=.025,right=.975,top=.94,bottom=.06)
ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
ax.text(.5,.97,'Verified rational 2-isogeny over F_(13^10)',ha='center',va='top',fontsize=12,weight='bold')
for x,key,color in [(.015,'source','#e8eff8'),(.575,'target','#e5f4ed')]:
    ax.add_patch(FancyBboxPatch((x,.30),.41,.51,boxstyle='round,pad=0.008',facecolor=color,edgecolor='#7589a8',lw=1))
    ax.text(x+.018,.755,'Source: Legendre short model' if key=='source' else 'Full rational 4-torsion target',fontsize=9,weight='bold',va='top')
    full=pair[key]['icv1']
    ax.text(x+.018,.65,textwrap.fill(full,width=41,break_long_words=True,break_on_hyphens=False),fontfamily='DejaVu Sans Mono',fontsize=7.6,va='top',linespacing=1.35)
ax.annotate('',xy=(.56,.55),xytext=(.445,.55),arrowprops={'arrowstyle':'->','color':'#235b9d','lw':2})
ax.text(.5,.65,'degree 2',ha='center',fontsize=9)
ax.text(.5,.39,'original-model\nkernel (0,0)',ha='center',fontsize=8)
ax.text(.5,.15,'Full ICV1 IDs name the short models; coordinate changes and map are recorded.',ha='center',fontsize=8)
ax.text(.5,.055,'Independent point counts agree; explicit point and target-halving replays pass.',ha='center',fontsize=8)
fig.savefig(ROOT/'verified_isogeny.svg');fig.savefig(ROOT/'verified_isogeny.pdf');plt.close(fig)
print('rendered two editable SVG figures and two companion PDFs')

# Keep generated vector files free of trailing whitespace.
for vector in (ROOT/'field_expansion.svg', ROOT/'verified_isogeny.svg'):
    vector.write_text('\n'.join(line.rstrip() for line in vector.read_text().splitlines())+'\n')
