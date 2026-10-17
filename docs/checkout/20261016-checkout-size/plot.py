# Lays out the recorded values in data.json; computes nothing new.
import json, pathlib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = pathlib.Path(__file__).parent
data = json.loads((here / "data.json").read_text())
rows = data["recipes"]
names = [r["label"] for r in rows]
tree = [r["worktree_bytes"] / 1e9 for r in rows]
git = [(r["git_bytes"] or 0) / 1e9 for r in rows]

fig, ax = plt.subplots(figsize=(7.5, 3.4))
y = range(len(rows))
ax.barh(y, tree, color="#4269d0", label="working tree")
ax.barh(y, git, left=tree, color="#a3b4e8", label=".git")
for i, r in enumerate(rows):
    total = tree[i] + git[i]
    note = f'{total:.2f} GB, {r["files"]:,} files'
    if r["clone_seconds"] is not None:
        note += f', {r["clone_seconds"]} s'
    ax.text(total + 0.08, i, note, va="center", fontsize=8)
ax.set_yticks(list(y), names, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("bytes on disk after clone (GB, 1e9)")
ax.set_xlim(0, max(t + g for t, g in zip(tree, git)) * 1.45)
ax.legend(fontsize=8, loc="lower right", frameon=False)
ax.set_title(f'Checkout size by recipe, main at {data["commit"]} (n=1 clone each)', fontsize=9)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(here / "checkout-size.svg")
fig.savefig(here / "checkout-size.png", dpi=150)
