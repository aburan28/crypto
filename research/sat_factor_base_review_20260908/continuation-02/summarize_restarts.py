"""Summarize the retained deterministic two-orbit perturbation runs."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
summary = []
for path in sorted(ROOT.glob("full-orbit-search-*-kick-*.jsonl")):
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    start = next(row for row in rows if row.get("kind") == "full_orbit_search_start")
    result = next(row for row in rows if row.get("kind") == "full_orbit_search_result")
    group = next(row for row in rows if row.get("kind") == "affine_quotient_base")
    summary.append(
        {
            "artifact": path.name,
            "objective": result["objective"],
            "kick_seed": start["kick_seed"],
            "initial_covered_target_orbits": start["covered_target_orbits"],
            "accepted_exchanges": result["accepted_exchanges"],
            "objective_covered_target_orbits": result["covered_target_orbits"],
            "covered_target_orbits_at_most_three": group["hits_at_most_three"],
            "covered_target_orbits_exactly_three": group["hits_exactly_three"],
            "selected_indices": result["selected_indices"],
        }
    )

summary.sort(key=lambda row: (row["objective"], row["kick_seed"]))
exact = [row for row in summary if row["objective"] == "exactly_three"]
hybrid = [row for row in summary if row["objective"] == "at_most_three"]
assert len(exact) == 9 and len(hybrid) == 17
assert max(row["objective_covered_target_orbits"] for row in exact) == 6256
assert max(row["objective_covered_target_orbits"] for row in hybrid) == 6300
assert [
    row["kick_seed"]
    for row in exact
    if row["objective_covered_target_orbits"] == 6256
] == [9118]
assert [
    row["kick_seed"]
    for row in hybrid
    if row["objective_covered_target_orbits"] == 6300
] == [9118]

output = {
    "kind": "full_orbit_restart_summary",
    "candidate_orbits": 6909,
    "selected_orbits": 4,
    "two_orbit_kick": True,
    "exactly_three_runs": len(exact),
    "at_most_three_runs": len(hybrid),
    "best_exactly_three": 6256,
    "best_at_most_three": 6300,
    "best_seed": 9118,
    "runs": summary,
    "scope": "deterministic restart census; not an exhaustive global search over four-orbit subsets",
}
(ROOT / "restart-summary.json").write_text(json.dumps(output, indent=2) + "\n")
print(json.dumps({key: output[key] for key in output if key != "runs"}))
