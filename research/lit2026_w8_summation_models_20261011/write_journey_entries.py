"""Write the W8 journey entries (layout only: copies recorded values, including the
ratios the Rust harness computed in compare.json).

usage: python -I write_journey_entries.py <repo root> <pr ref>
"""
import json
import sys
from pathlib import Path

repo = Path(sys.argv[1])
pr = sys.argv[2]
study = repo / "research/lit2026_w8_summation_models_20261011"
out = repo / "research/lit2026_speedup_program/journey/entries"
out.mkdir(parents=True, exist_ok=True)
cmp = json.loads((study / "reports/run-20261011T2131Z/compare.json").read_text())["result"]["cells"]
cells = {
    "K1/F2^15 m=2": ("k1-n15", "icv1-f2m15-t275-b7f03703", "r = 211"),
    "K0/F2^15 m=2": ("k0-n15", "icv1-f2m15-tm275-2d22ff5d", "r = 751"),
    "K1/F2^17 m=2": ("k1-n17", "icv1-f2m17-tm101-00378d4e", "r = 65587"),
    "K0/F2^23 m=2": ("k0-n23", "icv1-f2m23-t5197-69e76b73", "r = 2095853"),
    "K1/F2^23 m=2": ("k1-n23", "icv1-f2m23-tm5197-1f85e9e1", "r = 4196903"),
}
base = {"schema": "lit2026-journey-entry/v1", "workstream": "W8", "date": "2026-10-11", "pr": pr}
DET = "none (deterministic counter, identical in all 5 rounds)"
for cell, (tag, slug, sub) in cells.items():
    c = cmp[cell][0]
    arms = {a["arm"]: a for a in c["arms"]}
    rat = {(r["numerator"], r["denominator"]): r for r in c["ratios"]}
    ro = rat[("orig/sym", "cand/eg")]
    rp = rat[("cand/sym", "cand/eg")]
    e = dict(base)
    e.update({
        "slug": f"eg-quotient-s3-{tag}",
        "papers": ["ePrint 2026/1299"],
        "kind": "speedup",
        "workload": {
            "id": f"koblitz-fghr-sym-s3-m2-{tag}",
            "description": f"descended FGHR-symmetrised S3 (m = 2) decomposition solve, all roots, MatrixF4 cap 4, "
                           f"8 targets in the prime subgroup ({sub}) of {slug}",
            "curve_ids": [slug],
        },
        "metric": "F4 word XORs per round (sum over 8 target solves)",
        "unit": "counted ops",
        "isolation": "L0",
        "samples": 5,
        "arms": [
            {"role": "original", "label": "repo symmetrised system (sym), harness on the 4b6e09b78 library",
             "commit": "4b6e09b78", "value": arms["orig/sym"]["word_ops_totals_seen"][0]},
            {"role": "previous", "label": "sym in the candidate binary", "commit": "lit2026/w8-summation-models",
             "value": arms["cand/sym"]["word_ops_totals_seen"][0]},
            {"role": "candidate", "label": "S3 of E/<T> in w = 1/X (eg): sym with the parity bit eliminated",
             "commit": "lit2026/w8-summation-models", "value": arms["cand/eg"]["word_ops_totals_seen"][0]},
        ],
        "ratio_original_over_candidate": {"value": round(ro["word_ops_ratio"], 4),
                                          "interval": [round(ro["word_ops_ratio"], 4)] * 2, "interval_kind": DET},
        "ratio_previous_over_candidate": {"value": round(rp["word_ops_ratio"], 4),
                                          "interval": [round(rp["word_ops_ratio"], 4)] * 2, "interval_kind": DET},
        "outputs_identical": bool(ro["root_digests_identical"] and rp["root_digests_identical"]),
        "run_report": "research/lit2026_w8_summation_models_20261011/reports/run-20261011T2131Z/report.pdf",
        "note": (f"Wall time (L0, load 100-341, ratio of per-round totals): original/candidate median "
                 f"{ro['median_ratio_of_median_ms']:.2f} [{ro['per_round_ratio_min']:.2f}-{ro['per_round_ratio_max']:.2f}], "
                 f"previous/candidate {rp['median_ratio_of_median_ms']:.2f} "
                 f"[{rp['per_round_ratio_min']:.2f}-{rp['per_round_ratio_max']:.2f}]. Toy decomposition stage only; "
                 f"lifting an eg root tries both parities; no IC-vs-rho claim."),
    })
    (out / f"20261011-W8-eg-quotient-s3-{tag}.json").write_text(json.dumps(e, indent=2) + "\n")

R1 = "research/lit2026_w8_summation_models_20261011/reports/run-20261011T2017Z/report.pdf"
diag = [
    ("eg-quotient-s4-m3", "diagnostic", ["ePrint 2026/1299"], "koblitz-fghr-sym-s4-m3",
     "descended FGHR-symmetrised S4 (m = 3) against S4 of E/<T>, u-plain and x-direct on icv1-f2m15-t275-b7f03703 "
     "(r = 211), icv1-f2m15-tm275-2d22ff5d (r = 751), icv1-f2m9-t5-81e744be (r = 127)",
     ["icv1-f2m15-t275-b7f03703", "icv1-f2m15-tm275-2d22ff5d", "icv1-f2m9-t5-81e744be"], 40,
     "Res_s(sym, s^2+s+Sum w) = S4 of E/<T> in w (certificate IDC1hcf83176ae9a93715). At n = 15 the E/<T> system "
     "(Boolean degree 6, 12 vars) needs 1.55-1.89x fewer word ops than sym (degree 4, 13 vars) but 4.9-7.5x more "
     "splits and 1.6-2.0x the wall time; identical root sets. K1/F2^17 m = 3 stopped unfinished after 2 h."),
    ("wu-state-model-s3", "diagnostic", ["arXiv 2608.01675"], "binary-s3-coordinate-models",
     "descended S3 (m = 2) in Weierstrass x, Wu Kummer ratio k = x/d, Wu u = d/(x+d), binary Edwards w, on "
     "icv1-f2m15-tm275-2d22ff5d, icv1-f2m15-t5-d9fb183e, icv1-f2m17-t101-e6c4b64d, icv1-f2m17-t209-0eba1a1e "
     "(uniform affine targets)",
     ["icv1-f2m15-tm275-2d22ff5d", "icv1-f2m15-t5-d9fb183e", "icv1-f2m17-t101-e6c4b64d", "icv1-f2m17-t209-0eba1a1e"], 40,
     "Wu u is a Mobius frame of x and equals the FGHR frame (u -> u+1 is translation by T); k is x rescaled. S3 "
     "terms/total degree: x and k 5/4, u 11/6. Word ops u vs x: 1.2-1.4x fewer at n = 15, 3.1x more at n = 17; "
     "binary Edwards w is an E/<T> coordinate, 5.7-7.5x more. No model uniformly cheaper; root sets equal brute "
     "force everywhere."),
    ("frobenius-aux-variable", "negative", ["ePrint 2026/1927"], "koblitz-direct-x-s3-s4-frobenius-aux",
     "direct descended S_{m+1} on the Frobenius-stable x-base, with y = x^(2^k) auxiliary blocks, "
     "icv1-f2m15-t275-b7f03703, icv1-f2m15-tm275-2d22ff5d, icv1-f2m9-t5-81e744be",
     ["icv1-f2m15-t275-b7f03703", "icv1-f2m15-tm275-2d22ff5d", "icv1-f2m9-t5-81e744be"], 40,
     "Over the F2 descent x -> x^(2^k) is linear: same Boolean degree, same ideal after linear elimination. "
     "MatrixF4: 57-58x (m = 2), 701x (m = 3) more word ops; InheritedF4 (linear elimination): 1.03x (m = 2), "
     "8.6x (m = 3). Root sets identical."),
    ("convolution-e-splitting", "negative", ["arXiv 2609.18866"], "binary-s4-e-variable-split",
     "symmetrised S4 (b = 1) in e-variables with correspondence, V = <1..z^(l-1)>, n = 3l, l = 3, 4, 5, "
     "random target abscissae (polynomial-level decomposition, no curve group)",
     [], 40,
     "Convolution split (f = X1 X2): degree 3 -> 2, 4-14% fewer monomials, same split tree, 1.21-2.07x more word "
     "ops and 2.4-2.9x wall time. Lucas factorisation has no counterpart for F_{2^n}-valued summands; m = 4 "
     "measured at correspondence level only (no binary symmetrised S5 in the repo)."),
]
for slug, kind, pp, wid, desc, ids, n, note in diag:
    e = dict(base)
    e.update({
        "slug": slug, "papers": pp, "kind": kind,
        "workload": {"id": wid, "description": desc, "curve_ids": ids},
        "metric": "F4 word XORs per solve (median over targets x rounds)", "unit": "counted ops",
        "isolation": "L0", "samples": n, "arms": [],
        "ratio_original_over_candidate": None, "ratio_previous_over_candidate": None,
        "outputs_identical": True, "run_report": R1, "note": note,
    })
    (out / f"20261011-W8-{slug}.json").write_text(json.dumps(e, indent=2) + "\n")
print(sorted(p.name for p in out.glob("20261011-W8-*.json")))
