"""Summarise the 2000-permutation replays of cells that crossed a band in the
main run (PROTOCOL.md: POSITIVE needs p < 0.05/31 on a fresh replay).
Writes replay_summary.json.  Usage: python3 -I .../analysis/replay.py
"""
import glob
import json
import os

import numpy as np

ITEM = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
KEYS = ["maxZ", "Ztmax", "M1", "M10", "M100", "ZtM1", "ZtM10", "ZtM100"]
ALPHA = 0.05 / 31


def main():
    out = {}
    cells = sorted({os.path.basename(f).split("_replay")[0] for f in glob.glob(os.path.join(ITEM, "raw", "*_replay1*.json"))})
    for c in cells:
        vals = {k: [] for k in KEYS}
        obs = None
        for f in sorted(glob.glob(os.path.join(ITEM, "raw", c + "_replay1*.json"))):
            d = json.load(open(f))
            obs = d["observed"]
            for k in KEYS:
                vals[k] += d["perm_values"][k]
        row = {}
        for k in KEYS:
            v = np.asarray(vals[k])
            hi = (1 + (v >= obs[k]).sum()) / (1 + len(v))
            lo = (1 + (v <= obs[k]).sum()) / (1 + len(v))
            row[k] = dict(obs=obs[k], n=len(v), p_upper=float(hi), p_lower=float(lo),
                          q001=float(np.quantile(v, 0.001)), q999=float(np.quantile(v, 0.999)),
                          survives=bool(min(hi, lo) < ALPHA))
        out[c] = row
        print(c, {k: (round(row[k]["p_upper"], 4), round(row[k]["p_lower"], 4), row[k]["survives"]) for k in KEYS})
    json.dump(dict(alpha=ALPHA, cells=out), open(os.path.join(ITEM, "replay_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
