#!/usr/bin/env python3
"""Regenerate `docs/ic/calibration.json`: the unit's pinned conversion ratios.

usage: boundary_pin_calibration.py <run.json> [<run.json> …] [--out <path>] [--dry]

The boundary ledger's unit is group-addition equivalents, so every
non-addition count is divided by `ns_per_add` after being multiplied by
its own `ns_per_<unit>`.  Through Round 3 both factors were measured on
the host at the start of each run, and the *ratio* between them drifted:
on one machine, across three ladders, a median of 1.08 and up to 3.70 for
the same instance and unit.  That repriced rows which had done identical
native work by up to eight per cent, and it put a floor under every
cross-run comparison the ledger makes — 111 of the 166 rows in the
Round-3 comparison moved without an operation changing.

`AGENTS.md` §6 keeps operation counts as the metric "because they survive
hardware".  A conversion re-measured per run does not survive it, so the
ratios live here instead, and `ic boundary` prices with these rather than
with what it measured.  It still measures, and still reports what it
measured, because that is the wall-clock practicality note and because a
host that no longer resembles this one should be visible.

This script takes the median ratio per (instance, unit) over the runs
given.  Prefer several runs on one host; the median is what makes a
single unlucky measurement not become the unit.  Run it again when the
ladder gains an instance, or when the reference host changes — and when
it does, say so in the note, because re-pinning reprices every row.
"""
import json, statistics, sys, collections

UNITS = [
    "ns_per_double", "ns_per_sqrt", "ns_per_as_solve", "ns_per_s4_pair",
    "ns_per_lookup", "ns_per_row_op", "ns_per_word_xor", "ns_per_legendre",
    "ns_per_inversion", "ns_per_frobenius", "ns_per_canon",
]

paths = [a for a in sys.argv[1:] if not a.startswith("--")]
out_path = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--out=")), "docs/ic/calibration.json")
dry = "--dry" in sys.argv
if not paths:
    sys.exit(__doc__)

acc = collections.defaultdict(lambda: collections.defaultdict(list))
hosts = set()
for path in paths:
    run = json.load(open(path))
    hosts.add(run["host"]["cpu"])
    for inst in run["ledger"]["instances"]:
        # The measured factors, whether or not this run priced with them.
        cal = inst.get("calibration_measured") or inst["calibration"]
        add = cal.get("ns_per_add", 0.0)
        if add <= 0:
            continue
        key = f"{inst['regime']}/{inst['curve']['name']}"
        for unit in UNITS:
            value = cal.get(unit)
            if value:
                acc[key][unit].append(value / add)

instances, spreads = collections.OrderedDict(), []
for key in sorted(acc):
    row = collections.OrderedDict()
    for unit in UNITS:
        values = acc[key].get(unit)
        if not values:
            continue
        row[unit] = round(statistics.median(values), 6)
        if len(values) >= 2:
            spreads.append((max(values) / min(values), key, unit))
    instances[key] = row

doc = collections.OrderedDict([
    ("what", "Conversion ratios of the boundary ledger's unit: nanoseconds per native "
             "operation divided by nanoseconds per group addition, per instance per regime. "
             "`ic boundary` prices with these instead of with the factors it measures on the "
             "host, so that two runs price identical native counts identically."),
    ("why", "Measured per run, the ratios drifted by a median of 1.08 and up to 3.70 between "
            "runs on one host, repricing rows that had done identical work by up to eight per "
            "cent. Operation counts survive hardware; a conversion re-measured per run does not."),
    ("unit", "ns_per_<native unit> / ns_per_add, dimensionless"),
    ("method", f"median over {len(paths)} frozen run(s)"),
    ("host", sorted(hosts)[0] if len(hosts) == 1 else sorted(hosts)),
    ("sources", [p.split("/")[-1] for p in paths]),
    ("fallback", "An instance with no entry here keeps the factors measured on the host for "
                 "that run, and the report names the units that fell back, so a row priced the "
                 "old way says so."),
    ("regenerate", "python3 docs/ic/tools/boundary_pin_calibration.py docs/ic/runs/<run>.json …"),
    ("instances", instances),
])

if dry:
    print(json.dumps(doc, indent=1)[:2000])
else:
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
        fh.write("\n")
    print(f"written {out_path}")

print(f"instances {len(instances)}, units {sum(len(v) for v in instances.values())}")
if spreads:
    values = [s for s, _, _ in spreads]
    print(f"ratios collapsed: median spread {statistics.median(values):.2f}, max {max(values):.2f}")
    print("widest five:")
    for spread, key, unit in sorted(spreads, reverse=True)[:5]:
        print(f"   {spread:6.2f}  {key}  {unit}")
