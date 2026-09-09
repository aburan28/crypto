"""Independent scalar-orbit read-back of the full-universe result."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
R = 262543
N = 19


def rows(name):
    return [json.loads(line) for line in (ROOT / name).read_text().splitlines()]


roots = [
    value
    for value in range(1, R)
    if (value * value - value + 2) % R == 0 and pow(value, N, R) == 1
]
assert len(roots) == 1
LAMBDA = roots[0]


def orbit(anchor):
    values = set()
    current = anchor
    for _ in range(N):
        values.add(current)
        values.add((-current) % R)
        current = current * LAMBDA % R
    assert current == anchor and len(values) == 2 * N
    return sorted(values)


scalar_orbit = [-1] * R
scalar_orbit[0] = 0
target_representatives = []
for scalar in range(1, R):
    if scalar_orbit[scalar] >= 0:
        continue
    index = len(target_representatives)
    target_representatives.append(scalar)
    for value in orbit(scalar):
        assert scalar_orbit[value] in (-1, index)
        scalar_orbit[value] = index
assert len(target_representatives) == 6909
assert all(index >= 0 for index in scalar_orbit)

hybrid_rows = rows("full-orbit-search-hybrid-kick-9118.jsonl")
exact_rows = rows("full-orbit-search-exact3-kick-9118.jsonl")
hybrid = next(row for row in hybrid_rows if row.get("kind") == "full_orbit_search_result")
exact = next(row for row in exact_rows if row.get("kind") == "full_orbit_search_result")
hybrid_group = next(row for row in hybrid_rows if row.get("kind") == "affine_quotient_base")
exact_group = next(row for row in exact_rows if row.get("kind") == "affine_quotient_base")

expected_indices = [200, 1913, 2481, 5643]
expected_scalars = [203, 2143, 2901, 9853]
expected_points = [[16795, 144921], [1315, 90425], [8461, 38471], [6685, 369649]]
for result in [hybrid, exact]:
    assert result["candidate_orbits"] == len(target_representatives)
    assert result["selected_indices"] == expected_indices
    assert result["selected_scalar_representatives"] == expected_scalars
    assert result["selected_point_representatives"] == expected_points
for index, scalar in zip(expected_indices, expected_scalars):
    assert target_representatives[index] == scalar

selected_orbits = [orbit(scalar) for scalar in expected_scalars]
exact_support = set()
for a in range(4):
    anchor = expected_scalars[a]
    for b in range(a, 4):
        for c in range(b, 4):
            for second in selected_orbits[b]:
                for third in selected_orbits[c]:
                    total = (anchor + second + third) % R
                    if total:
                        exact_support.add(scalar_orbit[total])

hybrid_support = set(exact_support)
for a in range(4):
    hybrid_support.add(expected_indices[a])
    anchor = expected_scalars[a]
    for b in range(a, 4):
        for second in selected_orbits[b]:
            total = (anchor + second) % R
            if total:
                hybrid_support.add(scalar_orbit[total])

assert len(exact_support) == 6256
assert len(hybrid_support) == 6300
assert hybrid_support - exact_support
assert len(hybrid_support - exact_support) == 44


def words(support):
    result = [0] * ((len(target_representatives) + 63) // 64)
    for index in support:
        result[index // 64] |= 1 << (index % 64)
    return result


assert hybrid["selected_support_words"] == words(hybrid_support)
assert exact["selected_support_words"] == words(exact_support)
assert hybrid["covered_target_orbits"] == len(hybrid_support)
assert exact["covered_target_orbits"] == len(exact_support)
for group in [hybrid_group, exact_group]:
    assert group["hits_at_most_three"] == len(hybrid_support)
    assert group["hits_exactly_three"] == len(exact_support)
    assert group["target_orbit_size"] == 38
    assert group["quotient_points"] == 152
    assert group["projected_signed_orbits"] == 4
    assert sum(group["representative_multiplicities"]) * 38 + group[
        "zero_sum_triples_including_identity"
    ] == group["symmetric_cube_total"]

for run, objective, score in [
    (hybrid_rows, "at_most_three", 6300),
    (exact_rows, "exactly_three", 6256),
]:
    local = next(row for row in run if row.get("kind") == "full_orbit_one_exchange_local_optimum")
    final_round = local["round"]
    scans = [
        row
        for row in run
        if row.get("kind") == "full_orbit_removal_scan" and row["round"] == final_round
    ]
    assert len(scans) == 4
    assert {row["remove_position"] for row in scans} == {0, 1, 2, 3}
    assert all(row["objective"] == objective for row in scans)
    assert all(row["candidates_scored"] == 6906 for row in scans)
    assert all(row["best_covered_target_orbits"] == score for row in scans)
    assert local["best_neighbor_covered_target_orbits"] == score
    assert local["proper_neighbors_scored"] == 27620
    assert local["evaluations_including_current_reinsertions"] == 27624

sat_log = (ROOT / "selected-base-sat-round-trip.log").read_text()
sat_diagnostics = (ROOT / "selected-base-sat-round-trip-time.txt").read_text()
assert "test result: ok. 1 passed; 0 failed" in sat_log
assert "solver_calls: 1, models: 1" in sat_diagnostics
assert "exhausted: false, spurious: 0" in sat_diagnostics
assert "conflicts: 5472" in sat_diagnostics

messages = [
    "Reconstructed all 6909 signed-Frobenius target orbits from the public characteristic-polynomial root.",
    "Independent modular enumeration verified 6256 exactly-three and 6300 at-most-three target orbits for the selected four-orbit base.",
    "Both Rust support bitsets match the independent sets bit for bit; the group-law oracle matches the counts and symmetric-cube total.",
    "The recorded terminal scans cover all four removals and report no improving one-orbit neighbor for either objective.",
    "The selected explicit domain passes the recorded native-XOR S4 round trip in one model with 5472 conflicts and no invalid or exhausted result.",
]
(ROOT / "verification.txt").write_text("\n".join(messages) + "\n")
print("\n".join(messages))
