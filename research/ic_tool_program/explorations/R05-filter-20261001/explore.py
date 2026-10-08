#!/usr/bin/env python3
"""The exploration before R05's declaration (README.md): the prototype's
presence filter against the build it sits on, at three suite sizes. It
is not a round's measurement.

    IC_BASE=<ic at eaa842aa> IC_CAND=<ic at 0318f016> python3 explore.py run <out dir>
    tar -xJf runs.tar.xz && python3 explore.py summarise runs > summary.json

`run` interleaves the two arms on `M1`'s two rows at each size, three
rounds, ABAB, each process isolated by the harness. `summarise` reads
the run tree only.
"""
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "harness"))
import bench  # noqa: E402
import stats  # noqa: E402

SIZES = [(0, 41), (0, 53), (0, 61)]
ROUNDS = 3


def rows() -> list[dict]:
    out = [r for r in bench.slug_rows(bench.suite_rows("S"))
           if r["recipe_seed"] == 201 and (r["a"], r["n"]) in SIZES]
    assert len(out) == 6, len(out)
    return out


def run(out: Path) -> None:
    arms = {"base": Path(os.environ["IC_BASE"]), "cand": Path(os.environ["IC_CAND"])}
    bench.interleave(arms, rows(), ROUNDS, out)


def summarise(out: Path) -> dict:
    counts = {"processes": 0, "contended": 0, "failed": 0}
    for rec in out.rglob("*.isolation.jsonl"):
        for line in rec.read_text().splitlines():
            r = json.loads(line)["run"]
            counts["processes"] += 1
            counts["contended"] += bool(r["contended"])
            counts["failed"] += r["exit_status"] != 0
    sizes = {}
    for (a, n) in SIZES:
        cold, collect, same = [], [], True
        for r in (r for r in rows() if (r["a"], r["n"]) == (a, n)):
            for k in range(1, ROUNDS + 1):
                pa = bench.figure_path(out / "base" / r["id"] / f"r{k}.price.json")
                pb = bench.figure_path(out / "cand" / r["id"] / f"r{k}.price.json")
                if not (bench.clean(pa) and bench.clean(pb)):
                    continue
                ra, rb = bench.load(pa), bench.load(pb)
                same &= bench.outputs(ra) == bench.outputs(rb)
                cold.append(stats.ic_cold_ns(ra) / stats.ic_cold_ns(rb))
                collect.append(stats.setup_phase_ns(ra)["collect"] / stats.setup_phase_ns(rb)["collect"])
        logs = [math.log(x) for x in cold]
        mean = sum(logs) / len(logs)
        sd = math.sqrt(sum((x - mean) ** 2 for x in logs) / (len(logs) - 1))
        sizes[bench.curve_slug(a, n)] = {
            "cold_base_over_cand": stats.geo_ci(cold),
            "cold_pairs": [round(x, 4) for x in cold],
            "cold_log_sd": round(sd, 4),
            "collect_base_over_cand": stats.geo_ci(collect),
            "outputs_identical": same,
        }
    return {"what_this_is": "the exploration before R05's declaration, not a round's measurement: the "
                            "prototype filter (0318f016) over its base (eaa842aa), cold time and collection, "
                            "base over candidate",
            "accounting": counts, "sizes": sizes}


if __name__ == "__main__":
    step, out = sys.argv[1], Path(sys.argv[2]).resolve()
    if step == "run":
        run(out)
    elif step == "summarise":
        print(json.dumps(summarise(out), indent=1))
    else:
        raise SystemExit(__doc__)
