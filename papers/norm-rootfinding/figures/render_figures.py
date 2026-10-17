"""Presentation only: plot existing JSON counts; no algebra, oracle, or timing."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
HERE = Path(__file__).resolve().parent
ARMS = json.loads((HERE.parent/"data/E2a.json").read_text())["evidence"]["arms"]
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                     "svg.hashsalt":"rootfinding-paper-v1"})
def save(fig, name):
    for ext in ("pdf","svg","png"):
        fig.savefig(HERE/f"{name}.{ext}",bbox_inches="tight",dpi=150)
    plt.close(fig)
fig, ax = plt.subplots(figsize=(8,3.9))
for i, a in enumerate(ARMS):
    ax.bar(i-.18,a["external_variables"],.34,color="#49758b",
           label="External (measured)" if i==0 else None)
    val = a.get("native_internal_variables_after_flush",a["preclean_internal_estimate"])
    ax.bar(i+.18,val,.34,color="#cf8c3a" if a["native_constructed"] else "white",
           edgecolor="#9c631e",hatch=None if a["native_constructed"] else "///",
           label=("Internal (measured)" if i==0 else "Internal estimate; not loaded" if i==5 else None))
ax.set_yscale("log")
ax.set_ylim(10000,1500000)
ax.set_ylabel("Coordinates (log scale)")
ax.set_xticks(range(len(ARMS)),[str(a["support_limit"]) for a in ARMS])
ax.set_xlabel("Affine support limit")
ax.set_title("Native XOR cuts change the cost of affine contraction",loc="left",weight="bold",fontsize=12)
ax.grid(axis="y",alpha=.2)
ax.set_axisbelow(True)
ax.legend(fontsize=8,loc="upper left",frameon=False)
fig.text(.1,-.055,"Source: E2a, frozen native profile. No SAT solve; no measured speedup.",fontsize=9)
save(fig,"native_coordinates")
fig, ax = plt.subplots(figsize=(8.4,4))
ax.set_xlim(0,10);ax.set_ylim(0,5);ax.axis("off")
def box(x,y,w,h,text,color="#edf3f5",edge="#44616c",dash=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.08",
        facecolor=color,edgecolor=edge,linestyle="--" if dash else "-",linewidth=1.2))
    ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=9)
def arrow(a,b):
    ax.annotate("",xy=b,xytext=a,arrowprops={"arrowstyle":"->","color":"#44616c","lw":1.3})
box(.2,3.7,2.4,.9,"Finite field + curve\nfixed target + source domain")
box(3.25,3.7,2.9,.9,"Exact binary norm chart\nsources + r + kernel bit\nTheorem 2")
box(6.8,3.7,2.9,.9,"Exact affine compiler\nnative CDCL + XOR\nTheorem 3")
arrow((2.7,4.15),(3.15,4.15));arrow((6.25,4.15),(6.7,4.15))
box(6.8,1.95,2.9,.9,"SAT: replay literal norm\nblock full primary model\nUNSAT: chart exhaustion")
arrow((8.25,3.6),(8.25,2.95))
box(3.25,1.95,2.9,.9,"Budget / cap / UNKNOWN\nlower bound only\nexact counts unset",color="#fff0da",edge="#9c631e")
arrow((6.7,2.4),(6.25,2.4))
box(.2,.15,5.95,.95,"OPEN bridge: original S3 chain roots\nauxiliary fibers + ghosts + exceptional branches",color="#fff7ee",edge="#9c631e",dash=True)
arrow((6.8,1.95),(6.25,.85))
box(6.8,.15,2.9,.95,"Separate component:\nsupplied-H root extraction\nbinary / odd fields")
ax.text(.2,3.15,"Implemented chart contracts do not imply original-system completeness.",
        fontsize=9,color="#6c4c23")
save(fig,"architecture")
