"""Replay every source line and count its exact three-sum pencil support."""

import collections
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_adaptive_private_columns import add, cell, mul, neg, on_curve  # noqa: E402


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


paths = {
    "old_quotient": HERE / "type_i_line_guided_collector.json",
    "new_quotient": HERE / "type_i_line_full_new_quotient_points.json",
    "old_quartic": HERE / "type_i_line_base_points.tsv",
    "new_quartic": HERE / "type_i_line_new_points.tsv",
    "source_pairs": HERE / "type_i_line_new_point_pairs.tsv",
    "density": HERE / "type_i_character_density.json",
}
old = json.loads(paths["old_quotient"].read_text())
p = int(old["field_prime"])
old_q = [tuple(map(int, point)) for point in old["base_points_fp"]]
new_q = [tuple(map(int, point)) for point in
         json.loads(paths["new_quotient"].read_text())]
old_c = [tuple(map(int, line.split()[1:])) for line in
         paths["old_quartic"].read_text().splitlines()]
new_c = [tuple(map(int, line.split())) for line in
         paths["new_quartic"].read_text().splitlines()]
pairs = [tuple(map(int, line.split())) for line in
         paths["source_pairs"].read_text().splitlines()]
assert len(old_q) == len(old_c) == 141
assert len(new_q) == len(new_c) == len(pairs) == 9784
assert len(set(old_c + new_c)) == 9925
assert len(set(old_q + new_q)) == 9925
assert len(set(old_q + new_q) |
           {(x, (-y) % p) for x, y in old_q + new_q}) == 19850
assert all(on_curve(point, p) for point in old_q + new_q)


def collinear(a, b, c):
    return ((b[0] - a[0]) * (c[1] - a[1]) -
            (b[1] - a[1]) * (c[0] - a[0])) % p == 0


def dsu_find(parent, index):
    while parent[index] != index:
        parent[index] = parent[parent[index]]
        index = parent[index]
    return index


parent = list(range(141))
degrees = [0] * 141
edges = set()
line_endpoints = []
incident_lines = [[] for _ in range(141)]
line_section_target = None
star_targets = {}
for edge_index in range(4892):
    i, j = pairs[2 * edge_index]
    assert pairs[2 * edge_index + 1] == (i, j)
    assert 1 <= i < j <= 141
    assert (i, j) not in edges
    edges.add((i, j))
    line_endpoints.append((i - 1, j - 1))
    incident_lines[i - 1].append(edge_index)
    incident_lines[j - 1].append(edge_index)
    x, y = new_c[2 * edge_index:2 * edge_index + 2]
    assert x != y
    assert collinear(old_c[i - 1], old_c[j - 1], x)
    assert collinear(old_c[i - 1], old_c[j - 1], y)
    qx, qy = new_q[2 * edge_index:2 * edge_index + 2]
    section = add(add(add(old_q[i - 1], old_q[j - 1], p), qx, p), qy, p)
    if line_section_target is None:
        line_section_target = section
    assert section == line_section_target
    for center, other in ((i - 1, j - 1), (j - 1, i - 1)):
        target = add(add(old_q[other], qx, p), qy, p)
        if center in star_targets:
            assert star_targets[center] == target
        else:
            star_targets[center] = target
        degrees[center] += 1
    parent[dsu_find(parent, i - 1)] = dsu_find(parent, j - 1)

components = collections.Counter(dsu_find(parent, i) for i in range(141))
assert len(star_targets) == 141 and min(degrees) > 0
assert len(set(star_targets.values())) == 141
assert len(components) == 1
assert sum(degrees) == 9784
assert sum(degree - 1 for degree in degrees) == 9784 - 141

# A tree in the source-line adjacency graph gives an explicit independent
# collision basis: line e minus its earlier parent has a private new-point
# coordinate at e, so the rows are triangular in discovery order.
queue = collections.deque([0])
seen_lines = {0}
collision_tree = []
while queue:
    line = queue.popleft()
    for vertex in line_endpoints[line]:
        for neighbor in incident_lines[vertex]:
            if neighbor not in seen_lines:
                seen_lines.add(neighbor)
                collision_tree.append([neighbor, line, vertex + 1])
                queue.append(neighbor)
assert len(seen_lines) == 4892 and len(collision_tree) == 4891
previous = {0}
for child, ancestor, vertex in collision_tree:
    assert child not in previous and ancestor in previous
    assert vertex - 1 in line_endpoints[child]
    assert vertex - 1 in line_endpoints[ancestor]
    previous.add(child)
assert len(previous) == 4892
density = json.loads(paths["density"].read_text())
coefficients = [int(value) for value in density["sextic_coefficients_decimal"]]
r = int(old["subgroup_order"])
inv2 = pow(2, -1, r)
star_cell_by_vertex = {
    vertex: cell(neg(mul(inv2, target, p), p), coefficients, p)[::-1]
    for vertex, target in star_targets.items()
}
star_cells = collections.Counter(star_cell_by_vertex.values())
weighted_star_cells = collections.Counter()
for vertex, category in star_cell_by_vertex.items():
    weighted_star_cells[category] += degrees[vertex]
assert sum(star_cells.values()) == 141
assert sum(weighted_star_cells.values()) == 9784

# One private residual-point pair per source line makes the 4892 line-section
# vectors independent.  The connected line-adjacency graph gives 4891
# independent differences, and each star difference is one such difference.
out = {
    "status": "verified_full_source_pencil_graph",
    "field_prime": str(p),
    "old_classes": len(old_q),
    "new_classes": len(new_q),
    "source_lines": len(edges),
    "graph_components": len(components),
    "minimum_vertex_degree": min(degrees),
    "maximum_vertex_degree": max(degrees),
    "source_graph_degree_histogram": dict(sorted(collections.Counter(degrees).items())),
    "three_sum_decompositions": sum(degrees),
    "distinct_pencil_targets": len(star_targets),
    "pencil_targets_by_Q_H_cell": dict(sorted(star_cells.items())),
    "pencil_decompositions_by_Q_H_cell": dict(sorted(weighted_star_cells.items())),
    "duplicate_pencil_decompositions": sum(degrees) - len(star_targets),
    "independent_pencil_collision_rows": len(collision_tree),
    "collision_tree_rows_child_parent_shared_vertex": collision_tree,
    "line_section_target": [str(value) for value in line_section_target],
    "input_sha256": {name: sha(path) for name, path in paths.items()},
    "source_sha256": sha(Path(__file__)),
    "boundary": "These are exact source-line pencil collisions. Each line has a private new-point pair; these rows alone allow arbitrary old-column values.",
}
path = HERE / "type_i_line_full_pencil_graph.json"
path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
print("FULL_PENCIL_GRAPH_PASS", json.dumps({
    "lines": len(edges), "components": len(components),
    "three_sum_decompositions": out["three_sum_decompositions"],
    "distinct_targets": out["distinct_pencil_targets"],
    "collision_rank": out["independent_pencil_collision_rows"],
}, sort_keys=True))
