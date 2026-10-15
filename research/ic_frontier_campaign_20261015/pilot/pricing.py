"""Locate the count behind every unpriced counter of a record (recorded values only).

An unpriced name `X_uncharged` is counted in the record's top-level `counters`
or in a phase's `native` map, either as `X_uncharged` or as `X`.  Used by
summarise.py and the Phase 1 analysis; it never invents a count.
"""


def unpriced_counts(r):
    out = {}
    for name in r["cost"].get("unpriced") or []:
        base = name[: -len("_uncharged")] if name.endswith("_uncharged") else name
        if name in r.get("counters", {}):
            out[name] = r["counters"][name]
            continue
        total, found = 0, False
        for p in r.get("phases", []):
            for key in (name, base):
                if key in p.get("native", {}):
                    total += p["native"][key]
                    found = True
                    break
        if not found and base.startswith("solver_") and r.get("solver"):
            out[name] = r["solver"].get("ops")
            continue
        out[name] = total if found else None
    return out
