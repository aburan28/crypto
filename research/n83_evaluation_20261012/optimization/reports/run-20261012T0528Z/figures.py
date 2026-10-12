#!/usr/bin/env python3
"""Lay out recorded native values only; no research-result computation."""
from pathlib import Path
import json
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/n83-optimization-report-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"figures"
DATA=json.loads((ROOT/"data/statistics.json").read_text())
POLICY=json.loads((ROOT/"data/policy_selection.json").read_text())
NAVY,TEAL,BLUE,ORANGE,MUTED,GRID="#142B45","#087F83","#527AAB","#B66A32","#5B6A7A","#DAE4EB"


def save(fig,name):
    fig.savefig(OUT/f"{name}.svg",bbox_inches="tight")
    fig.savefig(OUT/f"{name}.png",dpi=210,bbox_inches="tight")
    plt.close(fig)


def policy_scores():
    rows=POLICY["candidate_scores"]
    fig,ax=plt.subplots(figsize=(10.0,6.1))
    fig.subplots_adjust(left=.24,right=.94,top=.82,bottom=.16)
    y=list(range(len(rows)-1,-1,-1))
    values=[r["primary_capped_arithmetic_mean_process_wall_s"] for r in rows]
    colors=[TEAL if r["candidate_id"]==POLICY["selected_candidate_id"] else BLUE for r in rows]
    ax.barh(y,values,color=colors,height=.65)
    ax.set_yticks(y,[r["candidate_id"] for r in rows])
    ax.set_xlim(0,5)
    ax.set_xlabel("Native arithmetic mean whole-process wall (s), 9 training cases per policy")
    ax.tick_params(axis="y",length=0)
    ax.grid(axis="x",color=GRID,lw=.7);ax.set_axisbelow(True)
    for pos,value in zip(y,values):
        ax.text(value+.06,pos,f"{value:.3f}",va="center",fontsize=9,color=MUTED)
    fig.text(.025,.96,"Registered Fold policy grid: complete training scores",fontsize=17.3,fontweight="bold",color=NAVY,va="top")
    fig.text(.025,.88,"All 12 policies complete and verify all 9 cases. Teal marks the native frozen selection.",fontsize=10.5,color=MUTED)
    fig.text(.025,.055,"3 training targets × 3 seeds; these are selection scores, not held-out estimates or a global optimum.",fontsize=9.6,color=MUTED)
    fig.text(.025,.012,"Source: native policy_selection.json · run-20261012T0528Z",fontsize=9.3,color=MUTED)
    save(fig,"policy_scores")


def holdout_rho():
    targets=DATA["per_target"]
    fig,ax=plt.subplots(figsize=(10,4.7));fig.subplots_adjust(left=.09,right=.97,top=.78,bottom=.25)
    x=list(range(5))
    for arm,shift,color,label in [("default",-.18,BLUE,"Default Wide/32/DP4/off"),("selected",.18,TEAL,"Selected Fold/64/DP6/off")]:
        ax.bar([a+shift for a in x],[t[arm]["mean_verified_actual_process_wall_s"] for t in targets],width=.32,color=color,label=label)
        for n,t in enumerate(targets):
            rows=[r for r in DATA["raw_rows"]["rho"] if r["target_id"]==t["target_id"] and r["arm"]==arm]
            offsets=[-.08,0,.08]
            ax.scatter([n+shift+offsets[j] for j in range(len(rows))],
                       [r["elapsed_process_wall_s"] for r in rows],color=NAVY,s=18,zorder=3,edgecolor="white",linewidth=.4)
    ax.set_xticks(x,[f"Q{n+1}" for n in x]);ax.set_ylim(bottom=0)
    ax.set_ylabel("Whole-process wall (s)");ax.set_xlabel("Held-out target block")
    ax.grid(axis="y",color=GRID,lw=.7);ax.set_axisbelow(True)
    ax.legend(frameon=False,fontsize=9,ncol=2,loc="upper left")
    fig.text(.025,.96,"Held-out rho: all 30 seed solves and native target means",fontsize=17,fontweight="bold",color=NAVY,va="top")
    fig.text(.025,.88,"Bars are native three-seed means; dots are all recorded individual process times.",fontsize=10.5,color=MUTED)
    fig.text(.025,.08,"Five target blocks under the orbit-exclusion law; no caps or failed answers. The policy was frozen before holdout.",fontsize=9.5,color=MUTED)
    fig.text(.025,.025,"Source: native statistics.json · run-20261012T0528Z",fontsize=9.3,color=MUTED)
    save(fig,"holdout_rho")


def ic_cases():
    targets=DATA["per_target"]
    fig,ax=plt.subplots(figsize=(10,4.7));fig.subplots_adjust(left=.10,right=.97,top=.78,bottom=.25)
    x=list(range(5))
    rho=[t["selected"]["mean_verified_online_ms"] for t in targets]
    ic=[t["prepared_ic"]["target_ms"] for t in targets]
    ax.bar([a-.17 for a in x],rho,width=.31,color=TEAL,label="Selected rho: native 3-seed mean")
    ax.bar([a+.17 for a in x],ic,width=.31,color=ORANGE,label="Prepared IC: one target interval")
    ax.set_xticks(x,[f"Q{n+1}" for n in x]);ax.set_ylabel("Complete prepared-online wall (ms)")
    ax.set_ylim(0,16000);ax.set_xlabel("Held-out target block")
    ax.grid(axis="y",color=GRID,lw=.7);ax.set_axisbelow(True)
    ax.legend(frameon=False,fontsize=9,ncol=2,loc="upper right")
    ax.text(.17,13500,"13,034.488 ms",ha="center",fontsize=9.5,color=ORANGE,fontweight="bold")
    fig.text(.025,.96,"All five IC cases: four short intervals and one tail",fontsize=17.3,fontweight="bold",color=NAVY,va="top")
    fig.text(.025,.88,"No case is removed. IC is prepared once; each rho target mean includes three fresh seed draws.",fontsize=10.4,color=MUTED)
    fig.text(.025,.08,"Geometric IC/rho ratio 0.4310 [0.1380,1.7764]; arithmetic mean-time ratio 1.5426 points the other way.",fontsize=9.7,color=MUTED)
    fig.text(.025,.025,"Native estimators; no definitive overall online winner. Source: statistics.json · run-20261012T0528Z",fontsize=9.4,color=MUTED)
    save(fig,"ic_cases")


def estimator_plot():
    fields=[("Default/selected rho, process",DATA["primary_default_over_selected_process_wall"],TEAL),
            ("Default/selected rho, online",DATA["secondary_default_over_selected_online"],BLUE),
            ("Prepared IC/selected rho, online",DATA["secondary_prepared_ic_over_selected_rho_online"],ORANGE)]
    fig,ax=plt.subplots(figsize=(10,4.3));fig.subplots_adjust(left=.32,right=.96,top=.76,bottom=.25)
    for y,(label,v,color) in zip([2,1,0],fields):
        center=v["ratio"];lo,hi=v["ci95_percentile"]
        # Errorbar endpoint differences are layout geometry, not new estimates.
        ax.errorbar(center,y,xerr=[[center-lo],[hi-center]],fmt="o",color=color,capsize=5,markersize=7,elinewidth=2)
    ax.set_yticks([2,1,0],[a[0] for a in fields]);ax.tick_params(axis="y",length=0)
    ax.set_xscale("log");ax.set_xlim(.09,2.6);ax.set_ylim(-.6,2.6)
    ticks=[.1,.2,.5,1,2];ax.xaxis.set_major_locator(FixedLocator(ticks));ax.xaxis.set_major_formatter(FixedFormatter([str(v) for v in ticks]));ax.minorticks_off()
    ax.axvline(1,color=NAVY,linestyle=(0,(4,3)),lw=1.3);ax.grid(axis="x",color=GRID,lw=.7)
    ax.set_xlabel("Native geometric mean target-block ratio (log axis)")
    fig.text(.025,.96,"Native held-out ratios and 95% target-block intervals",fontsize=16.8,fontweight="bold",color=NAVY,va="top")
    fig.text(.025,.87,"Five conditionally sampled targets; three rho seeds per arm and one prepared IC interval per target",fontsize=10.2,color=MUTED)
    fig.text(.025,.08,"Each label preserves its own numerator/denominator. IC's geometric interval includes parity; its arithmetic mean differs.",fontsize=9.2,color=MUTED)
    fig.text(.025,.025,"Source: native statistics.json · run-20261012T0528Z · no intervals estimated by this renderer",fontsize=9.2,color=MUTED)
    save(fig,"heldout_intervals")


if __name__=="__main__":
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.spines.top":False,"axes.spines.right":False,
                         "svg.fonttype":"none","savefig.facecolor":"white"})
    OUT.mkdir(exist_ok=True)
    policy_scores();holdout_rho();ic_cases();estimator_plot()
    print("Rendered native recorded scores/means/raw values/intervals only; no analysis or resampling.")
