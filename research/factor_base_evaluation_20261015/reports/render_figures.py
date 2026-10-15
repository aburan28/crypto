#!/usr/bin/env python3
"""Render recorded summary values only; does not run or analyze experiments."""
import json, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/factor-base-render-cache')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent
KINDS=['coset','low_x','low_digit_weight','random']
LABEL={'coset':'Coset','low_x':'Low x','low_digit_weight':'Low digit weight','random':'Random'}
def figures(summary):
    palette=['#1d6b73','#bf763b','#6675aa','#7f8b8c']
    plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'font.size':10})
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for j,kind in enumerate(KINDS):
        rows=[r for r in summary['generic'] if r['kind']==kind]
        x=np.arange(3)+(j-1.5)*.19
        axes[0].bar(x,[100*r['coverage'] for r in rows],width=.18,label=LABEL[kind],color=palette[j])
        axes[1].bar(x,[r['rank_per_cpu'] for r in rows],width=.18,color=palette[j])
    for ax in axes:ax.set_xticks(range(3),['16 bits','20 bits','24 bits']);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    axes[0].set_ylabel('Targets with a three-point decomposition (%)')
    axes[1].set_ylabel('New independent rows / total CPU second')
    axes[0].set_title('Exact coverage on common uniform targets');axes[1].set_title('Exact-join prototype: setup included')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False)
    fig.tight_layout(rect=(0,.09,1,1));fig.savefig(OUT/'base_comparison.svg');fig.savefig(OUT/'base_comparison.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,3.3))
    for i,kind in enumerate(KINDS[1:]):
        c=next(c for c in summary['contrasts'] if c['bits']==24 and c['control']==kind)['throughput']
        ax.errorbar(c['ratio'],i,xerr=[[c['ratio']-c['lower']],[c['upper']-c['ratio']]],fmt='o',capsize=4,color='#1d6b73')
    ax.axvline(1,color='#888888',linewidth=1);ax.set_yticks(range(3),['Coset / '+LABEL[k] for k in KINDS[1:]])
    ax.set_xlabel('Independent-row throughput ratio; 95% whole-field bootstrap interval')
    ax.set_title('24-bit paired comparisons (8 independent fields)');ax.grid(axis='x',alpha=.15)
    fig.tight_layout();fig.savefig(OUT/'paired_intervals.svg');fig.savefig(OUT/'paired_intervals.png',dpi=150);plt.close(fig)

if __name__ == '__main__':
    figures(json.loads((OUT.parent/'results/summary.json').read_text()))
