"""Copy recorded W12 values into lit2026 journey entries (layout only).

Usage: python -I write_journey.py <repo-root> <pr-ref>
Values are copied from the run's scaling.json and fits.json; nothing is
computed here.
"""
import json
import sys

root = sys.argv[1]
PR = sys.argv[2]
run = 'research/lit2026_w12_flagged_claims_20261011/reports/run-20261012T0009Z'
s = json.load(open(f'{root}/{run}/scaling.json'))
f = json.load(open(f'{root}/{run}/fits.json'))['fits_log2_vs_log2r']
pc = {(p['method'].startswith('zero'), p['bits']): p for p in s['per_curve']}
out = f'{root}/research/lit2026_speedup_program/journey/entries/'
for pr in s['paired_ratio']:
    b = pr['bits']
    z = pc[(True, b)]
    r = pc[(False, b)]
    g = pr['ratio_zero_minor_over_rho_geomean']
    lo, hi = pr['ci95_bootstrap']
    e = {"schema": "lit2026-journey-entry/v1", "workstream": "W12", "slug": f"zero-minor-vs-rho-fp{b}",
         "date": "2026-10-11", "pr": PR,
         "papers": ["arXiv 2607.09814", "arXiv 2310.04132", "arXiv 1703.07544"], "kind": "negative",
         "workload": {"id": f"zero-minor-planted-fp{b}",
                      "description": f"8 planted ECDLP targets on {z['slug']} (find_prime_order_curve({b}, 59297), prime order r = {z['r']}, cofactor 1, whole group)",
                      "curve_ids": [z['slug']]},
         "metric": "median counted field ops per solved target (F_p mul+sqr+inv)", "unit": "counted ops",
         "isolation": "L0", "samples": z['runs'],
         "arms": [{"role": "original", "label": "program baseline (method absent)", "commit": "4b6e09b78", "value": None},
                  {"role": "reference", "label": "rho.negation (rho_walk_with + NegationClasses), 32 runs", "commit": "c70d57ec7", "value": r['ops_median']},
                  {"role": "candidate", "label": z['method'] + " (arXiv 2607.09814 as specified)", "commit": "c70d57ec7", "value": z['ops_median']}],
         "ratio_original_over_candidate": {"value": None, "interval": None, "interval_kind": "none"},
         "ratio_zero_minor_over_reference": {"value": g, "interval": [lo, hi],
                                             "interval_kind": "paired bootstrap 95% (geometric mean over targets of zero-minor/rho field ops)"},
         "outputs_identical": True,
         "run_report": run + "/report.pdf",
         "note": f"All 8 zero-minor and 32 rho runs recovered the planted logarithm; the zero-minor method costs {g:.0f}x rho.negation per target in counted field ops (95% bootstrap {lo:.0f}-{hi:.0f}). A replay reproduced every count."}
    json.dump(e, open(out + f"20261011-W12-zero-minor-vs-rho-fp{b}.json", "w"), indent=2)


def fl(k):
    v = f[k]
    return {"value": v['slope'], "interval": v['slope_ci95_bootstrap'], "n": v['n']}


e = {"schema": "lit2026-journey-entry/v1", "workstream": "W12", "slug": "zero-minor-slope", "date": "2026-10-11", "pr": PR,
     "papers": ["arXiv 2607.09814"], "kind": "diagnostic",
     "workload": {"id": "zero-minor-ladder-fp16-fp26",
                  "description": "48 planted targets, 8 on each registered prime-order ladder curve find_prime_order_curve(b, 59297), b = 16..26 step 2, whole group",
                  "curve_ids": [p['slug'] for p in s['curves']]},
     "metric": "OLS slope of log2(counted field ops) against log2 r", "unit": "counted ops", "isolation": "L0",
     "samples": 48, "arms": [],
     "ratio_original_over_candidate": {"value": None, "interval": None, "interval_kind": "none"},
     "slopes_bootstrap95": {"rho.negation field ops": fl("rho.negation [field ops]"),
                            "zero-minor all phases": fl("zero_minor, d by Eq.(4) rule [field ops, all phases]"),
                            "zero-minor search phases": fl("zero_minor, d by Eq.(4) rule [field ops, search phases only]"),
                            "zero-minor d-subsets covered": fl("zero_minor, d by Eq.(4) rule [d-subsets covered]")},
     "outputs_identical": True, "run_report": run + "/report.pdf",
     "note": "Falsifier met: zero-minor search-phase slope 0.905 (95% bootstrap 0.769-1.021) >= 0.5, and a constant 260x-6715x rho; candidate minors examined ~ r (median covered/r 0.37-2.1 per curve). No sub-rho behaviour, consistent with the paper's own Table 2."}
json.dump(e, open(out + "20261011-W12-zero-minor-slope.json", "w"), indent=2)
e = {"schema": "lit2026-journey-entry/v1", "workstream": "W12", "slug": "hu-eprint-2025-015", "date": "2026-10-11", "pr": PR,
     "papers": ["ePrint 2025/015", "ePrint 2024/1923"], "kind": "negative",
     "workload": {"id": "hu-iica-not-implementable",
                  "description": "no workload: the full text of ePrint 2025/015 was not obtainable (Cloudflare challenge, not bypassed); the abstract leaves IICA undefined and concerns DLP in F_p^*, not the ECDLP",
                  "curve_ids": []},
     "metric": "none", "unit": "counted ops", "isolation": "L0", "samples": 0, "arms": [],
     "ratio_original_over_candidate": {"value": None, "interval": None, "interval_kind": "none"},
     "outputs_identical": None, "run_report": "research/lit2026_w12_flagged_claims_20261011/REPORT.pdf",
     "note": "Blocked by access, not by mathematics: the algorithm is unspecified in every reachable source. The companion ePrint 2024/1923 (abstract) states no subexponential prime-field ECDLP method is known."}
json.dump(e, open(out + "20261011-W12-hu-eprint-2025-015.json", "w"), indent=2)
