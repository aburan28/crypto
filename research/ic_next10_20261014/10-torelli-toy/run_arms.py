"""Complete arms: hyperelliptic g = 2..8, plane quartics, trigonal quintics,
smooth plane quintics.  Writes data/<arm>.json = {models, smooth, hist}."""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fam_hyper import enumerate_hyper
from fam_plane import enumerate_plane

HERE = os.path.dirname(os.path.abspath(__file__))


def save(name, hist, smooth, models, secs, extra=None):
    d = {"arm": name, "models_enumerated": models, "smooth_models": smooth,
         "seconds": round(secs, 1),
         "hist": {",".join(map(str, k)): v for k, v in sorted(hist.items())}}
    if extra:
        d.update(extra)
    json.dump(d, open(os.path.join(HERE, "data", name + ".json"), "w"), indent=0)
    print(name, models, smooth, len(hist), f"{secs:.0f}s", flush=True)


def main(which):
    for arm in which:
        t = time.time()
        if arm.startswith("hyp"):
            g = int(arm[3:])
            h, tot = enumerate_hyper(g)
            save(arm, h, tot, 2 ** ((g + 2) + (2 * g + 3)) - 2 ** (2 * g + 3), time.time() - t,
                 {"genus": g, "group_order": 6 * 2 ** (g + 2)})
        elif arm == "quartic":
            h, ns, nt = enumerate_plane(4, 3)
            save(arm, h, ns, nt, time.time() - t, {"genus": 3, "group_order": 168})
        elif arm == "trigonal5":
            h, ns, nt = enumerate_plane(5, 5, trigonal=True)
            save(arm, h, ns, nt, time.time() - t, {"genus": 5, "group_order": 24})
        elif arm == "quintic":
            h, ns, nt = enumerate_plane(5, 6)
            save(arm, h, ns, nt, time.time() - t, {"genus": 6, "group_order": 168})


if __name__ == "__main__":
    main(sys.argv[1:])
