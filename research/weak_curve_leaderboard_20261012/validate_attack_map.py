#!/usr/bin/env python3
"""Check WAM1 attack-index integrity and the deployed static data copy."""

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SITE = HERE.parents[1] / "website" / "weakcurves" / "dist"


def main() -> None:
    criteria = json.loads((HERE / "criteria.json").read_text())
    attack = json.loads((HERE / "attack_map.json").read_text())
    assert attack["version"] == "WAM1"
    assert attack["assessment"] == "unassessed_on_indexed_models"
    assert (HERE / "attack_map.json").read_bytes() == (SITE / "attack_map.json").read_bytes()
    assert (HERE / "criteria.json").read_bytes() == (SITE / "criteria.json").read_bytes()
    audit = json.loads((HERE / "ghs_genus128_cost_audit.json").read_text())
    assert audit["schema"] == "weakcurves.ghs_genus128_cost_audit.v1"
    assert audit["producer_language"] == "rust"
    assert (HERE / "ghs_genus128_cost_audit.json").read_bytes() == \
        (SITE / "ghs_genus128_cost_audit.json").read_bytes()
    for source, published in (("snapshot.json", "snapshot.json"),
                              ("isogeny_graph.json", "isogeny_graph.json"),
                              ("IDENTIFIERS.md", "identifiers.md")):
        assert (HERE / source).read_bytes() == (SITE / published).read_bytes(), source

    ids = [entry["id"] for entry in criteria] + [entry["id"] for entry in attack["items"]]
    assert len(ids) == len(set(ids)), "attack and criterion IDs must be distinct"
    assert "ic_relation_support_rank" in ids, "relation-support criterion missing"
    assert "ic_scalar_anchor_overlap" in ids, "scalar-anchor criterion missing"
    assert len(attack["items"]) >= 30, "the expansion index lost entries"
    for item in criteria:
        for field in ("id", "name", "method", "domain", "recognition", "witness",
                      "construction", "status", "source", "url"):
            assert isinstance(item.get(field), str) and item[field].strip(), (item["id"], field)
        assert re.fullmatch(r"[a-z][a-z0-9_]*", item["id"]), item["id"]
        assert item["url"].startswith("https://"), item["url"]
    decomposition_ids = {
        "hyperelliptic_divisor_decomposition",
        "hyperelliptic_large_prime",
        "hyperelliptic_high_genus",
        "jacobian_isogeny_decomposition",
        "elliptic_quotient_prym",
        "jv_intermediate_field_cover",
        "jv_quadratic_sieve",
    }
    assert decomposition_ids <= set(ids), "hyperelliptic decomposition gate missing"
    for item in attack["items"]:
        for field in ("id", "name", "family", "scope", "requires", "test", "source", "url"):
            assert isinstance(item.get(field), str) and item[field].strip(), (item["id"], field)
        assert re.fullmatch(r"[a-z][a-z0-9_]*", item["id"]), item["id"]
        assert item["url"].startswith("https://"), item["url"]

    graph = json.loads((HERE / "isogeny_graph.json").read_text())
    nodes = {node["id"]: node for node in graph["nodes"]}
    for edge in graph["edges"]:
        assert edge["degree"] > 1
        assert edge["source"] in nodes and edge["target"] in nodes
        assert nodes[edge["source"]]["class_id"] == nodes[edge["target"]]["class_id"]
        assert edge["source_ref"] and edge["evidence"].startswith("verified_")
    assert graph["counts"]["graph_nodes"] == len(nodes)
    print(f"WAM1 OK: {len(criteria)} indexed criteria + {len(attack['items'])} research checks; "
          f"{len(nodes)} graph nodes, {len(graph['edges'])} certified directed maps")


if __name__ == "__main__":
    main()
