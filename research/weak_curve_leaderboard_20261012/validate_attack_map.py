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

    ids = [entry["id"] for entry in criteria] + [entry["id"] for entry in attack["items"]]
    assert len(ids) == len(set(ids)), "attack and criterion IDs must be distinct"
    assert len(attack["items"]) >= 30, "the expansion index lost entries"
    decomposition_ids = {
        "hyperelliptic_divisor_decomposition",
        "hyperelliptic_large_prime",
        "hyperelliptic_high_genus",
        "jacobian_isogeny_decomposition",
        "elliptic_quotient_prym",
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
