#!/usr/bin/env python3
"""Summarize the matched-arithmetic recheck's raw JSON outputs.

Reads the three producer JSONL outputs (original KS, matched-arith KS, IC)
plus the two callgrind instruction-count logs, and prints the numbers this
note's table needs: wall seconds, field_ops totals (KS arms only), callgrind
retired-instruction totals, and the ratios against the pre-registered rule.
Nothing here recomputes a number the producers/callgrind did not themselves
report; this only extracts and divides.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load_summary(path, kind):
    summary = None
    fixtures = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get("kind") == kind:
                summary = row
            else:
                fixtures.append(row)
    return summary, fixtures


def callgrind_total_instructions(annotate_path):
    """Parse `callgrind_annotate`'s totals line, e.g.
    '1,234,567,890  PROGRAM TOTALS' at the top of the summary block."""
    text = Path(annotate_path).read_text()
    m = re.search(r"([\d,]+)\s+PROGRAM TOTALS", text)
    if not m:
        return None
    return int(m.group(1).replace(",", ""))


def main():
    ks_orig_path = HERE / "corpus_n53-ks-growing-1024-v1.jsonl"
    ks_matched_path = HERE / "matched_arith_n53_L1024.jsonl"
    ic_path = HERE / "ic_n53_K440.jsonl"

    if ks_orig_path.exists():
        s, fx = load_summary(ks_orig_path, "rho_ks_batch_summary")
        print("== original KS (unmodified), n=53 L=1024 ==")
        if s:
            print("  in_process_ms:", s.get("in_process_ms"))
            print("  total_walk_steps:", s.get("total_walk_steps"))
            print("  all_verified:", s.get("all_verified"))
        print("  fixtures parsed:", len(fx))

    if ks_matched_path.exists():
        s, fx = load_summary(ks_matched_path, "rho_ks_batch_summary")
        print("== matched-arithmetic KS, n=53 L=1024 ==")
        if s:
            print("  in_process_ms:", s.get("in_process_ms"))
            print("  total_walk_steps:", s.get("total_walk_steps"))
            print("  all_verified:", s.get("all_verified"))
            fo = s.get("field_ops", {})
            print("  field_ops.total_field_ops:", fo.get("total_field_ops"))
            print("  field_ops.total_mul_calls:", fo.get("total_mul_calls"))
            print("  field_ops.total_sqr_calls:", fo.get("total_sqr_calls"))
            print("  field_ops.total_inv_calls:", fo.get("total_inv_calls"))
        print("  fixtures parsed:", len(fx))

    if ic_path.exists():
        s, fx = load_summary(ic_path, "compact_orbit_dlp_summary")
        print("== IC (koblitz_orbit_dlp_fast), n=53 K=440 ==")
        if s:
            print("  timing_ms:", json.dumps(s.get("timing_ms"), indent=2))
            print("  rank_attempts:", s.get("rank_attempts"))
            print("  rank_relations:", s.get("rank_relations"))
            print("  rank_failures:", s.get("rank_failures"))
            print("  targets_solved / targets:", s.get("targets_solved"), "/", s.get("targets"))
        print("  fixtures parsed:", len(fx))

    for name, path in [
        ("matched-arith KS", HERE / "callgrind_ks_matched.annotate.txt"),
        ("IC", HERE / "callgrind_ic.annotate.txt"),
    ]:
        if path.exists():
            total = callgrind_total_instructions(path)
            print(f"== callgrind total instructions, {name}: {total}")


if __name__ == "__main__":
    main()
