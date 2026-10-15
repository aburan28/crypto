#!/usr/bin/env python3
"""Bind hyperelliptic route properties to exact WeakCurves identities."""

import hashlib
import json
import math
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "decomposition_properties.json"
SITE_OUTPUT = ROOT / "website/weakcurves/dist/decomposition_properties.json"


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trim(poly):
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def remainder(a, b, p):
    a = trim(a[:])
    b = trim(b[:])
    if not b:
        raise ZeroDivisionError("polynomial division by zero")
    inverse = pow(b[-1], -1, p)
    while len(a) >= len(b):
        shift = len(a) - len(b)
        scale = a[-1] * inverse % p
        for i, coefficient in enumerate(b):
            a[i + shift] = (a[i + shift] - scale * coefficient) % p
        trim(a)
    return a


def squarefree(poly, p):
    a = trim(poly[:])
    b = trim([(i * a[i]) % p for i in range(1, len(a))])
    while b:
        a, b = b, remainder(a, b, p)
    return len(a) == 1


def multiply(a, b, p):
    result = [0] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            result[i + j] = (result[i + j] + left * right) % p
    return trim(result)


def formula_cover(record):
    """Replay the prime quadratic pullback from docs/curves/COVERS.md."""
    assert record["degree"] == 1 and int(record["p"]) > 3
    p, a = int(record["p"]), int(record.get("a", record.get("model_a")))
    b = int(record.get("b", record.get("model_b")))
    assert (4 * a**3 + 27 * b**2) % p != 0, record["id"]
    c = next(c for c in range(4) if (c**3 + a * c + b) % p)
    f = [(c**3 + a * c + b) % p, 0, (3 * c * c + a) % p,
         0, (3 * c) % p, 0, 1]
    assert squarefree(f, p), record["id"]
    x = [c, 0, 1]
    replay = multiply(multiply(x, x, p), x, p)
    replay[0] = (replay[0] + a * c + b) % p
    replay[2] = (replay[2] + a) % p
    assert replay == f, record["id"]
    hexes = lambda values: [hex(value) for value in values]
    return {
        "status": "formula_replayed",
        "field_relation": "same_field",
        "degree": 2,
        "genus": 2,
        "certificate": {"construction": "prime_quadratic_pullback_v1",
                        "degree": 2, "genus": 2, "f": hexes(f), "h": [],
                        "x": hexes([c, 0, 1]), "y_0": [], "y_v": ["0x1"]},
        "field_check": "source_prime_assumption",
        "subfield_descent": "not_tested",
        "subgroup_transfer": "not_tested",
        "source": "docs/curves/COVERS.md; build_decomposition_properties.py",
    }


def branch(record, label_key, screen_key):
    if label_key not in record:
        return {"status": "unassessed", "label": None, "screen": None}
    label, screen = record[label_key], record[screen_key]
    if label == "yes":
        status = "positive_model_label"
    elif screen == "outside_recorded_domain":
        status = "outside_recorded_domain"
    elif screen in ("excluded_by_depth_condition", "excluded_by_model_conditions"):
        status = "excluded_by_screen"
    elif label == "unclassified_in_this_index":
        status = "necessary_screen_passed"
    else:
        status = "no_recorded_model"
    return {"status": status, "label": label, "screen": screen}


def ghs_route(route):
    invariants, destination = route["invariants"], route["destination"]
    return {"evidence": route["evidence"], "transport_screen": route["transport"],
            "base_degree": invariants.get("base_degree"),
            "relative_degree": invariants.get("relative_degree"),
            "magic_number": invariants.get("magic_number"),
            "magic_number_bound": invariants.get("magic_number_bound"),
            "type_i": invariants.get("type_i"),
            "relative_trace_b_zero": invariants.get("relative_trace_b_zero"),
            "genus_projection": destination.get("genus"),
            "destination_field_bits": destination.get("field_bits"),
            "formula_class": route.get("class"),
            "formula_log2_total": route.get("cost", {}).get("log2_total")}


DEFINITIONS = {
    "same_field_cover": {"label": "Same-field hyperelliptic cover", "test": "Exact H to E equation, genus, map degree and declared field; same-field geometry does not test subfield descent."},
    "cubic_cover": {"label": "Odd-field cubic genus-three branch", "test": "Exact model label and necessary conductor-depth screen in the recorded cubic domain."},
    "quadratic_cover": {"label": "Odd-field quadratic genus-three branch", "test": "Exact model label and necessary model conditions in the recorded quadratic domain."},
    "ghs_descent": {"label": "Binary GHS lower-field route", "test": "Magic number, projected genus, base-field degree and transport-size screen; a cover equation and subgroup map remain separate."},
    "jv_intermediate_field_cover": {"label": "Joux–Vitse lower-field cover", "test": "For F_(p^6), inspect both H/F_(p²) genus-three and H/F_(p³) genus-two branches; construct the actual lower-field H and test norm–pullback on the selected subgroup."},
    "jv_quadratic_sieve": {"label": "Joux–Vitse quadratic sieve", "test": "For the odd-characteristic H/F_(p²) branch, verify the bihomogeneous saturated ideal, root sieve, full and partial relation yield, rank and target descent."},
    "hyperelliptic_divisor_decomposition": {"label": "Divisor factor-base decomposition", "test": "Destination H, divisor encoding, frozen factor base, natural full-relation yield, rank, final linear algebra and target descent."},
    "hyperelliptic_large_prime": {"label": "Large-prime divisor relations", "test": "Partial-relation yield, verified recombination graph, rank, final linear algebra and target descent."},
    "hyperelliptic_high_genus": {"label": "High-genus smoothness route", "test": "Exact destination genus and field, factor-base degree bound, natural-query smoothness and complete cost."},
    "jacobian_isogeny_decomposition": {"label": "Jacobian factorization", "test": "Weil polynomial over the destination field, factorization, explicit projection and selected-subgroup image."},
    "elliptic_quotient_prym": {"label": "Elliptic quotient and complement", "test": "Cover geometry, complementary factor, computable Jacobian projection and kernel on the selected subgroup."},
    "subgroup_transport": {"label": "Selected-subgroup transport", "test": "Exact selected order, generator and target, map images, kernel intersection and scalar compatibility."},
    "route_cost": {"label": "Complete route cost", "test": "Cover search, map construction, divisor or factor solve, target descent and independent replay against the same-group reference."},
}


def build():
    snapshot_path = HERE / "snapshot.json"
    covers_path = ROOT / "docs/curves/covers.json"
    registry_path = ROOT / "docs/curves/registry.json"
    transfers_path = ROOT / "docs/curves/jacobian-transfers.json"
    snapshot, covers = read(snapshot_path), read(covers_path)
    registry, transfers = read(registry_path), read(transfers_path)
    # These reports were built at different registry revisions. Bind each
    # imported certificate by exact ICV1 identity and full model hash below.
    by_id = {item["icv1"]: item for item in registry["curves"]}
    by_slug = {item["slug"]: item for item in covers["curves"]}
    by_transfer = {item["icv1"]: item for item in transfers["curves"]}
    assert len(by_slug) == len(covers["curves"])
    rows = snapshot["prime_examples"] + snapshot["models"] + snapshot["named_curves"]
    assert len(rows) == len({row["id"] for row in rows})
    records = {}
    for row in rows:
        identity = row["id"]
        checked_cover = by_slug.get(row["slug"])
        if checked_cover:
            source_model = by_id[identity]
            model_hash = hashlib.sha256(source_model["model_json"].encode()).hexdigest()
            assert checked_cover["model_sha256"] == model_hash
            assert source_model["slug"] == row["slug"]
            if row.get("sha256"):
                assert row["sha256"] == model_hash
            assert checked_cover["status"] == "verified_over_declared_field"
            assert checked_cover["field_relation"] == "same_field"
            cert = checked_cover["certificate"]
            cover = {"status": "verified_geometry", "field_relation": "same_field",
                     "degree": cert["degree"], "genus": cert["genus"],
                     "certificate": cert, "model_sha256": model_hash,
                     "field_check": checked_cover["field_check"],
                     "subfield_descent": checked_cover["subfield_descent"],
                     "subgroup_transfer": checked_cover["subgroup_transfer"],
                     "source": "docs/curves/covers.json"}
        else:
            cover = formula_cover(row)
        transfer = by_transfer.get(identity)
        routes = [ghs_route(route) for route in transfer["routes"]
                  if route["kind"] == "ghs_weil_descent"] if transfer else []
        exact_routes = [route for route in routes if route["evidence"] == "exact_invariant"]
        ghs_status = ("exact_invariant" if exact_routes else
                      "heuristic_projection" if routes else "unassessed")
        high_genus = [route for route in exact_routes
                      if route["transport_screen"] == "size_admits"
                      and route["genus_projection"] is not None
                      and int(route["genus_projection"]) >= 4]
        same_field_formula = next((route for route in transfer["routes"]
                                   if route["kind"] == "same_field_cover"), None) if transfer else None
        selected_order = row.get("generator_order") or row.get("r")
        candidate_order = transfer["subgroup"]["order"] if transfer else None
        sextic_odd = row["degree"] == 6 and int(row["p"]) % 2 == 1
        cubic_label = branch(row, "cubic", "cubic_screen") if sextic_odd else None
        jv_tower_status = ("cubic_model_label" if cubic_label and
                           cubic_label["status"] == "positive_model_label" else
                           "tower_screen_only" if sextic_odd else "outside_screened_towers")
        records[identity] = {
            "field": {"characteristic": row["p"], "degree": row["degree"],
                      "bits": row["bits"]},
            "same_field_cover": cover,
            "cubic_cover": branch(row, "cubic", "cubic_screen"),
            "quadratic_cover": branch(row, "quadratic", "quadratic_screen"),
            "ghs_descent": {"status": ghs_status, "routes": routes,
                            "constructed_cover": None},
            "jv_intermediate_field_cover": {
                "status": jv_tower_status,
                "towers": ([{"intermediate_degree": 2, "relative_degree": 3,
                             "destination_genus": 3, "geometry": "hyperelliptic",
                             "model_status": cubic_label["status"],
                             "model_screen": cubic_label["screen"]},
                            {"intermediate_degree": 3, "relative_degree": 2,
                             "destination_genus": 2, "geometry": "hyperelliptic",
                             "model_status": "unassessed_scholten"}] if sextic_odd else []),
                "lower_field_cover": None, "norm_pullback_kernel": None,
                "selected_subgroup_image": None},
            "jv_quadratic_sieve": {
                "status": ("cubic_model_label" if jv_tower_status == "cubic_model_label"
                           else "tower_screen_only" if sextic_odd
                           else "outside_screened_towers"),
                "intermediate_degree": 2 if sextic_odd else None,
                "relation_length": 8 if sextic_odd else None,
                "saturated_ideal": None, "nontrivial_full_relations": None,
                "partial_relations": None, "independent_rank": None,
                "target_descent": None},
            "hyperelliptic_divisor_decomposition": {
                "status": "unmeasured", "factor_base_digest": None,
                "natural_queries": None, "full_relations": None,
                "independent_rank": None, "target_descent": None},
            "hyperelliptic_large_prime": {
                "status": "unmeasured", "partial_relations": None,
                "collision_graph": None, "recombined_rank": None},
            "hyperelliptic_high_genus": {
                "status": "genus_projection_only" if high_genus else "unassessed",
                "candidate_routes": high_genus, "smoothness_measurement": None},
            "jacobian_isogeny_decomposition": {
                "status": "weil_polynomial_open", "destination_field": "source_field",
                "weil_polynomial": None, "factorization": None,
                "explicit_projection": None},
            "elliptic_quotient_prym": {
                "status": "cover_geometry_only", "elliptic_factor_dimension": 1,
                "complement_dimension": cover["genus"] - 1,
                "explicit_jacobian_projection": None},
            "subgroup_transport": {
                "status": "unverified", "selected_order": str(selected_order) if selected_order else None,
                "gcd_degree_selected_order": math.gcd(cover["degree"], int(selected_order)) if selected_order else None,
                "registry_subgroup_order": candidate_order,
                "gcd_degree_registry_order": math.gcd(cover["degree"], int(candidate_order)) if candidate_order else None,
                "generator_image": None, "target_image": None,
                "kernel_intersection": None},
            "route_cost": {
                "status": "formula_only" if same_field_formula else "unmeasured",
                "same_field_formula_class": same_field_formula.get("class") if same_field_formula else None,
                "same_field_formula_log2_total": same_field_formula.get("cost", {}).get("log2_total") if same_field_formula else None,
                "measured_complete_cost": None},
        }
    counts = {key: dict(sorted(Counter(record[key]["status"] for record in records.values()).items()))
              for key in DEFINITIONS}
    assert len(records) == 922 and counts["same_field_cover"] == {
        "formula_replayed": 5, "verified_geometry": 917}
    return {"schema": "weakcurves.decomposition_properties.v1",
            "sources": {"snapshot_sha256": digest(snapshot_path),
                        "covers_sha256": digest(covers_path),
                        "registry_sha256": digest(registry_path),
                        "transfers_sha256": digest(transfers_path)},
            "definitions": DEFINITIONS, "counts": counts,
            "records": dict(sorted(records.items()))}


def main():
    result = json.dumps(build(), sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n"
    if sys.argv[1:] == ["--check"]:
        assert OUTPUT.read_text() == result, "research property export is stale"
        assert SITE_OUTPUT.read_text() == result, "site property export is stale"
    elif not sys.argv[1:]:
        OUTPUT.write_text(result)
        SITE_OUTPUT.write_text(result)
    else:
        raise SystemExit("usage: build_decomposition_properties.py [--check]")
    print("Decomposition properties OK: 922 exact identities; 917 checked cover certificates, 5 formula replays")


if __name__ == "__main__":
    main()
