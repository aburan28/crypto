#!/usr/bin/env python3
"""Independent consistency checks for the per-curve decomposition index."""

import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def load(path):
    return json.loads(path.read_text())


def main():
    path = HERE / "decomposition_properties.json"
    site = ROOT / "website/weakcurves/dist/decomposition_properties.json"
    assert path.read_bytes() == site.read_bytes()
    props = load(path)
    snapshot = load(HERE / "snapshot.json")
    covers = load(ROOT / "docs/curves/covers.json")
    registry = load(ROOT / "docs/curves/registry.json")
    transfer = load(ROOT / "docs/curves/jacobian-transfers.json")
    assert props["schema"] == "weakcurves.decomposition_properties.v1"
    for name, source in (("snapshot", HERE / "snapshot.json"),
                         ("covers", ROOT / "docs/curves/covers.json"),
                         ("registry", ROOT / "docs/curves/registry.json"),
                         ("transfers", ROOT / "docs/curves/jacobian-transfers.json")):
        assert props["sources"][name + "_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    rows = snapshot["models"] + snapshot["named_curves"] + snapshot["prime_examples"]
    by_slug = {item["slug"]: item for item in covers["curves"]}
    by_id = {item["icv1"]: item for item in registry["curves"]}
    transfers = {item["icv1"]: item for item in transfer["curves"]}
    assert len(props["records"]) == len(rows) == 922
    assert set(props["records"]) == {row["id"] for row in rows}
    for row in rows:
        record = props["records"][row["id"]]
        assert set(record) == {"field", *props["definitions"]}
        assert record["field"]["characteristic"] == row["p"]
        assert record["field"]["degree"] == row["degree"]
        cover = record["same_field_cover"]
        assert cover["field_relation"] == "same_field"
        assert cover["certificate"]["genus"] == cover["genus"]
        assert cover["certificate"]["degree"] == cover["degree"]
        if row["slug"] in by_slug:
            source_cover = by_slug[row["slug"]]
            source_curve = by_id[row["id"]]
            assert cover["status"] == "verified_geometry"
            assert source_cover["certificate"] == cover["certificate"]
            assert source_cover["model_sha256"] == hashlib.sha256(source_curve["model_json"].encode()).hexdigest()
            assert cover["subgroup_transfer"] == "not_tested"
        else:
            assert cover["status"] == "formula_replayed"
            assert row["degree"] == 1 and int(row["p"]) > 3
            assert cover["genus"] == cover["degree"] == 2
        if "cubic" in row:
            assert record["cubic_cover"]["label"] == row["cubic"]
            assert record["quadratic_cover"]["label"] == row["quadratic"]
            assert record["cubic_cover"]["screen"] == row["cubic_screen"]
            assert record["quadratic_cover"]["screen"] == row["quadratic_screen"]
        if row["id"] not in transfers:
            assert record["route_cost"]["status"] == "unmeasured"
        assert record["subgroup_transport"]["status"] == "unverified"
        assert record["jacobian_isogeny_decomposition"]["explicit_projection"] is None
        assert record["hyperelliptic_divisor_decomposition"]["independent_rank"] is None
        assert record["hyperelliptic_large_prime"]["recombined_rank"] is None
        assert record["route_cost"]["measured_complete_cost"] is None
    for key in props["definitions"]:
        expected = dict(sorted(Counter(row[key]["status"] for row in props["records"].values()).items()))
        assert props["counts"][key] == expected
    assert props["counts"]["cubic_cover"]["positive_model_label"] == 14
    assert props["counts"]["quadratic_cover"]["positive_model_label"] == 10
    assert props["counts"]["cubic_cover"]["necessary_screen_passed"] == 247
    assert props["counts"]["cubic_cover"]["excluded_by_screen"] == 209
    assert props["counts"]["ghs_descent"]["exact_invariant"] == 49
    for name, base in (("Oakley Group 3", 31), ("Oakley Group 4", 37)):
        row = next(row for row in snapshot["named_curves"] if row.get("name") == name)
        routes = props["records"][row["id"]]["ghs_descent"]["routes"]
        assert any(route["evidence"] == "exact_invariant" and route["base_degree"] == base
                   and route["relative_degree"] == 5 and route["magic_number"] == 5
                   and route["genus_projection"] == "15" for route in routes)
    print("Decomposition index OK: 922 identities, 917 source certificates, 5 formula replays, 11 property filters")


if __name__ == "__main__":
    main()
