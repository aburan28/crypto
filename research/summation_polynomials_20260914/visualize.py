#!/usr/bin/env python3
"""Create portable charts from frozen measurements; never rerun benchmarks."""
import csv
import os
from html import escape
import json
import hashlib
import shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
OUT=REPO/'docs/performance-gains'
OUT.mkdir(exist_ok=True)
VIZ=Path(os.environ['CRYPTO_VISUALIZATION_DIR']) if os.environ.get('CRYPTO_VISUALIZATION_DIR') else None
a=json.loads((HERE/'run-01/comparison.json').read_text())
b=json.loads((HERE/'confirmation-01/comparison.json').read_text())
assert a['complete'] and b['complete']
variants=['reference','setup','fixed','horner','combined']
labels=['Reference','Setup reuse only','Fixed slots','Grouped Horner','Combined · retained']
colors=['#93a1a4','#aa8b5b','#69879e','#637bad','#167f73']
local={(x['variant'],x['metric']):x for x in a['aggregates'] if x['kind']=='micro'}
full={(x['variant'],x['metric']):x for x in b['aggregates']}
ref=json.loads((HERE/'run-01/micro-p67-seed11-r0-reference.json').read_text())
new=json.loads((HERE/'run-01/micro-p67-seed11-r0-combined.json').read_text())
counts=[x['specialization_fp_muls']//(x['targets']*x['evaluation_repetitions']) for x in [ref,new]]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelcolor':'#34494f',
 'text.color':'#173139','xtick.color':'#53656b','ytick.color':'#173139','axes.edgecolor':'#d4dad9',
 'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'axes.titleweight':'bold',
 'svg.fonttype':'none','pdf.fonttype':42,'savefig.facecolor':'#f7f6f2'})
fig,axes=plt.subplots(2,2,figsize=(15.8,11),facecolor='#f7f6f2')
fig.subplots_adjust(left=.13,right=.95,top=.77,bottom=.20,wspace=.57,hspace=.73)
fig.text(.065,.945,'Where the performance gains actually are',fontsize=25,weight='bold')
fig.text(.065,.908,'ECC index calculus  /  Summation polynomials  /  Frozen measurements · 14 September 2026',fontsize=11,color='#607178')
fig.text(.065,.867,'55.4% less specialization time',fontsize=16,weight='bold',color='#167f73')
fig.text(.53,.867,'Full-DLP improvement remains unproven',fontsize=14,weight='bold',color='#8d6735')
for ax in axes.flat:
 ax.set_facecolor('#f7f6f2');ax.grid(axis='x',color='#e2e5df',lw=.7);ax.set_axisbelow(True);ax.tick_params(axis='y',length=0,pad=10)

def local_chart(ax,metric,title,sub,zoom=False):
 y=np.arange(len(variants));ratios=np.array([local[v,metric]['ratio']*100 for v in variants]);ci=np.array([local[v,metric]['paired_95'] for v in variants])*100
 if not zoom:ax.barh(y,ratios,height=.54,color=colors,zorder=2)
 for i,v in enumerate(variants):
  ax.errorbar(ratios[i],y[i],xerr=[[ratios[i]-ci[i,0]],[ci[i,1]-ratios[i]]],fmt='o' if zoom else 'none',
   color=colors[i] if zoom else '#18363b',markersize=6,capsize=4,elinewidth=1.5,zorder=4)
  x=(max(ratios[i],ci[i,1])+.1) if zoom else (max(ratios[i],ci[i,1])+2)
  ax.text(x,y[i],f'{ratios[i]:.2f}%' if zoom else f'{ratios[i]:.1f}%',va='center',fontsize=10,weight='bold' if v=='combined' else 'normal')
 ax.set_yticks(y,labels);ax.invert_yaxis();ax.axvline(100,color='#567079',ls='--',lw=1,zorder=3)
 ax.set_xlim((98.3,101.35) if zoom else (0,117))
 ax.set_title(title,loc='left',fontsize=13,pad=33)
 ax.text(0,1.065,sub,transform=ax.transAxes,fontsize=10,color='#63747a')
 ax.set_xlabel('Time relative to reference (%) · lower is better',labelpad=10,fontsize=10)
local_chart(axes[0,0],'specialization_s','A  ·  Per-target specialization','Reference = 100%; whiskers show paired 95% intervals')
local_chart(axes[0,1],'setup_s','B  ·  Once-per-curve symbolic setup','Zoomed scale highlights the small setup change',True)
ax=axes[1,0];y=np.arange(2)
ax.barh(y,counts,height=.48,color=[colors[0],colors[-1]])
ax.set_yticks(y,['Corrected reference','Grouped Horner']);ax.invert_yaxis();ax.set_xlim(0,1580)
for i,v in enumerate(counts):ax.text(v+30,i,f'{v:,}',va='center',fontsize=13,weight='bold')
ax.set_title('C  ·  Real arithmetic savings',loc='left',fontsize=13,pad=33)
ax.text(0,1.065,'Fp multiplications per specialization · p = 67, seed = 11',transform=ax.transAxes,fontsize=10,color='#63747a')
ax.set_xlabel('Counted multiplications · lower is better',labelpad=10,fontsize=10)
ax.text(0,-.45,f'{counts[0]-counts[1]} fewer products ({100*(1-counts[1]/counts[0]):.1f}%). Historical 1,770 estimate excluded:\nits correction was accounting, not an optimization gain.',transform=ax.transAxes,fontsize=9,color='#63747a',linespacing=1.5)
ax=axes[1,1]
for i,v in enumerate(['legacy','horner','combined']):
 for metric,dy,c,marker in [('wall_s',-.13,'#637bad','o'),('cpu_s',.13,'#167f73','D')]:
  d=full[v,metric];x=d['ratio']*100;lo,hi=np.array(d['paired_95'])*100
  ax.errorbar(x,i+dy,xerr=[[x-lo],[hi-x]],fmt=marker,color=c,capsize=4,ms=6,elinewidth=1.7,
   label=('Wall time' if metric=='wall_s' else 'CPU time') if i==0 else None)
ax.set_yticks(range(3),['Legacy control','Grouped Horner','Combined · retained']);ax.invert_yaxis();ax.set_xlim(90,106)
ax.axvline(100,color='#567079',ls='--',lw=1)
ax.set_title('D  ·  Full-DLP confirmation',loc='left',fontsize=13,pad=33)
ax.text(0,1.065,'All intervals cross 100%: no established gain or regression',transform=ax.transAxes,fontsize=10,color='#63747a')
ax.set_xlabel('Time relative to reference (%) · lower is better',labelpad=10,fontsize=10)
ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.02,-.43),ncol=2,fontsize=10)
fig.text(.065,.063,'Evidence: 342 ablation runs + 140 confirmation runs. Polynomial and scalar checks passed throughout.',fontsize=10,color='#53676c')
fig.text(.065,.039,'Paired 95% intervals cluster repetitions by input. Shared host; small corpus. Local gains are not additive or an asymptotic advance.',fontsize=9,color='#63747a')
for ext in ['png','svg','pdf']:
 fig.savefig(OUT/f'summary.{ext}',dpi=180)
 if ext == 'svg':
  svg=OUT/'summary.svg'
  svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
plt.close(fig)
history=json.loads((REPO/'research/gaudry_allocation_20260914/run-01/comparison.json').read_text())
gpu=json.loads((REPO/'gpu/macaulay/benchmarks/pipeline-20260914/run-01/summary.json').read_text())
data={'variants':variants,'labels':labels,'colors':colors,'local':a['aggregates'],'confirmation':b['aggregates'],
 'counts':counts,'history':history['aggregates'],'gpu':gpu['rows'],
 'localPairs':a['pairs'],'confirmationPairs':b['pairs'],
 'sourceHashes':{str(p.relative_to(REPO)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'run-01/comparison.json',HERE/'confirmation-01/comparison.json',REPO/'research/gaudry_allocation_20260914/run-01/comparison.json',REPO/'gpu/macaulay/benchmarks/pipeline-20260914/run-01/summary.json']}}
(OUT/'data.json').write_text(json.dumps(data,indent=2)+'\n')
template=(HERE/'visualization.html').read_text()
names=dict(zip(variants,labels)); names['legacy']='Legacy control'
def cell(row):
 return f"<td>{row['ratio']*100:.2f}%<small>95%: {row['paired_95'][0]*100:.2f}–{row['paired_95'][1]*100:.2f}%</small></td>"
local_rows=[]; full_rows=[]; export=[]
for v in ['legacy',*variants]:
 category='accounting control' if v=='legacy' else 'corrected reference' if v=='reference' else 'engineering'
 local_rows.append('<tr'+(' class="retained"' if v=='combined' else '')+'><th scope="row">'+escape(names[v])+'</th><td>'+category+'</td>'+cell(local[v,'specialization_s'])+cell(local[v,'setup_s'])+'</tr>')
 for metric in ['specialization_s','setup_s']:
  r=local[v,metric]; export.append(['local',v,metric,r['ratio']*100,r['paired_95'][0]*100,r['paired_95'][1]*100,'percent_of_corrected_reference'])
for v in ['reference','legacy','horner','combined']:
 full_rows.append('<tr'+(' class="retained"' if v=='combined' else '')+'><th scope="row">'+escape(names[v])+'</th>'+cell(full[v,'wall_s'])+cell(full[v,'cpu_s'])+'<td>'+('Reference' if v=='reference' else 'Includes no change')+'</td></tr>')
 for metric in ['wall_s','cpu_s']:
  r=full[v,metric]; export.append(['confirmation',v,metric,r['ratio']*100,r['paired_95'][0]*100,r['paired_95'][1]*100,'percent_of_corrected_reference'])
with (OUT/'comparisons.csv').open('w',newline='') as f:
 writer=csv.writer(f,lineterminator='\n'); writer.writerow(['experiment','variant','metric','relative_time_percent','paired_95_low_percent','paired_95_high_percent','unit']); writer.writerows(export)
html=template.replace('__FROZEN_DATA__',json.dumps(data).replace('</','<\\/')).replace('__LOCAL_TABLE__',''.join(local_rows)).replace('__FULL_TABLE__',''.join(full_rows))
(REPO/'docs/performance-gains.html').write_text(html)
if VIZ:
 VIZ.mkdir(parents=True,exist_ok=True)
 shutil.copy2(REPO/'docs/performance-gains.html',VIZ/'performance-gains.html')
 # Preserve relative assets so the HTML and its downloads work offline.
 shutil.copytree(OUT,VIZ/'performance-gains',dirs_exist_ok=True)
print('Created dashboard, PNG, SVG, PDF and frozen chart data')
