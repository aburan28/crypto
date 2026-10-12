#!/usr/bin/env python3
"""Build the evidence-graded WeakCurves graph from frozen local records.

Same (q, trace) gives an F_q-isogeny class by Tate's theorem. Only the
independently verified walk records below contribute explicit map edges.
"""

import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WALK = ROOT / "research/isogeny_walk_p192_20261004/runs/p192-ell61-radius2/isogeny_routes.json"
WALK_VERIFY = WALK.with_name("VERIFY.json")
LARGE = ROOT / "docs/curves/large-prime-walks-20261010.json"
LARGE_MAPS = ROOT / "research/p192_p224_hundred_large_prime_isogenies_20261010/evidence"
INPUT = HERE / "snapshot.json"
OUTPUT = HERE / "isogeny_graph.json"


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_check(checks, criterion, status, detail):
    current = checks.get(criterion)
    priority = {"condition": 4, "geometry": 3, "map": 3, "candidate": 2, "screen": 1}
    if current is None or priority[status] > priority[current["status"]]:
        checks[criterion] = {"status": status, "detail": detail}


def atlas_checks(record):
    checks = {}
    grades = {
        "exact_condition": "condition",
        "covering_model": "geometry",
        "exact_invariant_candidate": "candidate",
        "candidate": "candidate",
    }
    for finding in record.get("findings", []):
        add_check(checks, finding["criterion_id"], grades[finding["grade"]], finding["statement"])
    if record.get("embedding_degree") is not None:
        for criterion in ("mov_weil_pairing", "frey_ruck_tate_pairing"):
            add_check(checks, criterion, "candidate", "Exact selected-order embedding degree; pairing transport and destination cost remain open.")
    if record.get("cubic_screen") == "necessary_condition_passed":
        add_check(checks, "cubic_norm_one", "screen", "Necessary conductor-depth condition passed; an exact cover is still required.")
    if record.get("cubic") == "yes" or record.get("quadratic") == "yes":
        add_check(checks, "cover_and_decomposition", "candidate", "An exact covering-model label exists; subgroup transport and decomposition cost are open.")
    if record["degree"] == 1 and int(record["p"]) > 3:
        p = int(record["p"])
        a = record.get("model_a", record.get("a"))
        b = record.get("model_b", record.get("b"))
        if a is not None and b is not None:
            a = int(a, 0) if isinstance(a, str) else a
            b = int(b, 0) if isinstance(b, str) else b
            if a % p == 0 and b % p and p % 3 == 1:
                add_check(checks, "glv_endomorphism", "screen", "j=0 and p≡1 (mod 3): a nontrivial rational order-three automorphism exists; selected-subgroup action is unverified.")
            if b % p == 0 and a % p and p % 4 == 1:
                add_check(checks, "glv_endomorphism", "screen", "j=1728 and p≡1 (mod 4): a rational order-four automorphism exists; selected-subgroup action is unverified.")
    return checks


def build():
    snapshot, walk, large = read(INPUT), read(WALK), read(LARGE)
    receipt = read(WALK_VERIFY)
    assert receipt["result"] == "pass" and receipt["isogeny_routes_sha256"] == digest(WALK)
    assert receipt["curves_verified"] == len(walk["curve_nodes"])
    assert receipt["edges_verified"] == len(walk["edges"])
    assert large["status"] == "PASS" and large["new_distinct_maps"] == len(large["records"])
    nodes = {}
    by_slug = {}
    by_uid = {}
    for record in snapshot["prime_examples"] + snapshot["models"] + snapshot["named_curves"]:
        ident = record["id"]
        assert ident not in nodes
        node = {
            "id": ident,
            "label": record.get("name", record.get("slug")),
            "slug": record["slug"],
            "p": record["p"], "degree": record["degree"], "trace": str(record["trace"]),
            "bits": record["bits"], "atlas": True,
            "checks": atlas_checks(record),
        }
        nodes[ident] = node
        assert node["slug"] not in by_slug
        by_slug[node["slug"]] = ident

    def ensure_node(uid, slug, p, trace, bits, label=None):
        ident = by_slug.get(slug, uid)
        if ident not in nodes:
            nodes[ident] = {"id": ident, "label": label or slug, "slug": slug,
                            "p": str(p), "degree": 1, "trace": str(trace),
                            "bits": bits, "atlas": False, "checks": {}}
            by_slug[slug] = ident
        node = nodes[ident]
        assert node["p"] == str(p) and node["degree"] == 1 and node["trace"] == str(trace)
        assert uid not in by_uid or by_uid[uid] == ident
        by_uid[uid] = ident
        node["curve_uid"] = uid
        return ident

    p192 = next(n for n in nodes.values() if n["label"] == "P-192")
    p224 = next(n for n in nodes.values() if n["label"] == "P-224")
    for item in walk["curve_nodes"]:
        p = item["field_characteristic"]
        t = int(p) + 1 - int(item["group_order"])
        ident = ensure_node(item["curve_uid"], item["icv1_slug"], p, t, 192)
        nodes[ident]["walk_ref"] = item["ref"]
        nodes[ident]["walk_depth"] = item["depth"]
        nodes[ident]["volcano_levels"] = item["proved_volcano_levels"]
    assert p192["curve_uid"] == walk["curve_nodes"][0]["curve_uid"]

    for item in large["records"]:
        assert item["status"] == item["kernel_check"] == item["codomain_check"] == item["subgroup_check"] == item["exact_rational_map_check"] == "PASS"
        root = p192 if item["source_icv1"].startswith("icv1-fp192-") else p224
        prefix = "icv1-fp192-" if root is p192 else "icv1-fp224-"
        assert all(item[f"{side}_icv1"].startswith(prefix + "t" + root["trace"] + "-") for side in ("source", "target"))
        for prefix in ("source", "target"):
            ensure_node(item[f"{prefix}_curve_uid"], item[f"{prefix}_icv1"], root["p"], root["trace"], root["bits"])

    edges = []
    for item in walk["edges"]:
        assert item["status"] == "verified" and item["separable"]
        source, target = by_uid[item["source_curve_uid"]], by_uid[item["target_curve_uid"]]
        assert nodes[source]["p"] == nodes[target]["p"] and nodes[source]["trace"] == nodes[target]["trace"]
        edges.append({"id": item["id"], "source": source, "target": target,
                      "degree": item["degree"], "direction": item["direction"],
                      "evidence": "verified_kernel_velu", "source_ref": str(WALK.relative_to(ROOT))})
    step_by_curve = Counter()
    for item in large["records"]:
        source, target = by_uid[item["source_curve_uid"]], by_uid[item["target_curve_uid"]]
        assert nodes[source]["p"] == nodes[target]["p"] and nodes[source]["trace"] == nodes[target]["trace"]
        curve = "p192" if item["source_icv1"].startswith("icv1-fp192-") else "p224"
        step = step_by_curve[curve]
        step_by_curve[curve] += 1
        map_path = LARGE_MAPS / curve / f"step-{step:03d}" / curve / f"ell-{item['ell']}.json"
        assert digest(map_path) == item["certificate_sha256"], map_path
        edges.append({"id": f"large-{item['ell']}-{item['index']}-{item['certificate_sha256'][:12]}",
                      "source": source, "target": target, "degree": item["ell"],
                      "direction": "recorded", "evidence": "verified_rational_map",
                      "certificate_sha256": item["certificate_sha256"],
                      "map_ref": str(map_path.relative_to(ROOT)),
                      "source_ref": str(LARGE.relative_to(ROOT))})
    assert len(edges) == len({e["id"] for e in edges})

    members = defaultdict(list)
    for node in nodes.values():
        key = (node["p"], node["degree"], node["trace"])
        members[key].append(node["id"])
    classes = []
    for (p, degree, trace), ids in sorted(members.items(), key=lambda x: (int(x[0][0]), x[0][1], int(x[0][2]))):
        raw = f"{p}:{degree}:{trace}".encode()
        class_id = "FQI1:" + hashlib.sha256(raw).hexdigest()[:16]
        for ident in ids:
            nodes[ident]["class_id"] = class_id
            nodes[ident]["class_size"] = len(ids)
        classes.append({"id": class_id, "p": p, "degree": degree, "trace": trace,
                        "size": len(ids), "members": sorted(ids)})
    degrees = defaultdict(set)
    for edge in edges:
        degrees[edge["source"]].add(edge["degree"])
        degrees[edge["target"]].add(edge["degree"])
    for ident, ell in degrees.items():
        add_check(nodes[ident]["checks"], "isogeny_transfer", "map",
                  "Explicit certified isogeny map of degree(s) " + ", ".join(map(str, sorted(ell))) + "; destination weakness and route cost require separate assessment.")

    counts = Counter()
    for node in nodes.values():
        if node["atlas"]:
            for criterion in node["checks"]:
                counts[criterion] += 1
    return {
        "schema": "weakcurves-isogeny-graph-v1",
        "meaning": "Class membership follows the same finite-field size and Frobenius trace; it is not an explicit map. Solid edges are certified maps from the named source records.",
        "source_sha256": {str(path.relative_to(ROOT)): digest(path) for path in (INPUT, WALK, WALK_VERIFY, LARGE)},
        "counts": {"atlas_nodes": sum(n["atlas"] for n in nodes.values()),
                   "graph_nodes": len(nodes), "classes": len(classes),
                   "multi_model_classes": sum(c["size"] > 1 for c in classes),
                   "explicit_edges": len(edges), "atlas_checks": dict(sorted(counts.items()))},
        "nodes": sorted(nodes.values(), key=lambda n: n["id"]),
        "classes": classes,
        "edges": sorted(edges, key=lambda e: e["id"]),
    }


def main():
    data = build()
    encoded = json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"
    if len(sys.argv) == 2 and sys.argv[1] == "--check":
        assert OUTPUT.read_text() == encoded, "isogeny_graph.json is stale"
    elif len(sys.argv) == 1:
        OUTPUT.write_text(encoded)
    else:
        raise SystemExit("usage: build_isogeny_graph.py [--check]")
    print(json.dumps(data["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
