#!/usr/bin/env python3
"""Validate the non-authoritative ECDLP audit planning mock.

JSON Schema catches structure. This verifier catches semantic invariants that
must fail closed: the planning identity, exact eight-workstream population,
metrics, references, dependency acyclicity, source hashes, report artifacts,
and the rule that a mock cannot promote claims.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from jsonschema import Draft202012Validator, FormatChecker


REPO = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MOCK = HERE / "ECDLP_AUDIT_PLAN.mock.json"
SCHEMA = HERE / "ecdlp-audit-plan.schema.json"


def canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def expected_plan_uid(record: dict) -> str:
    preimage = {
        "schema_version": record["schema_version"],
        "title": record["title"],
        "as_of": record["as_of"],
        "workstream_ids": [workstream["id"] for workstream in record["workstreams"]],
    }
    digest = hashlib.sha256(canonical_bytes(preimage)).hexdigest()
    return f"urn:ecdlp-audit-plan:1:sha256:{digest}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def check_acyclic(edges: list[tuple[str, str]], errors: list[str]) -> None:
    graph: dict[str, set[str]] = {f"WS{i}": set() for i in range(1, 9)}
    for source, target in edges:
        graph[source].add(target)

    temporary: set[str] = set()
    permanent: set[str] = set()

    def visit(node: str, trail: list[str]) -> None:
        if node in permanent:
            return
        if node in temporary:
            errors.append("cyclic blocking dependency: " + " -> ".join(trail + [node]))
            return
        temporary.add(node)
        for neighbor in sorted(graph[node]):
            visit(neighbor, trail + [node])
        temporary.remove(node)
        permanent.add(node)

    for node in sorted(graph):
        visit(node, [])


def main() -> int:
    record = json.loads(MOCK.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    errors: list[str] = []

    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    for error in sorted(validator.iter_errors(record), key=lambda item: list(item.path)):
        location = ".".join(str(part) for part in error.absolute_path) or "$"
        errors.append(f"schema {location}: {error.message}")

    require(record["plan_uid"] == expected_plan_uid(record), "plan_uid does not match the canonical planning preimage", errors)
    require(record["record_state"] == "mock_non_authoritative", "record must remain a non-authoritative mock", errors)
    require(record["authoritative"] is False, "mock cannot be authoritative", errors)
    require(record["promotion_allowed"] is False, "mock cannot allow promotion", errors)

    workstreams = record["workstreams"]
    workstream_ids = [item["id"] for item in workstreams]
    expected_ids = [f"WS{i}" for i in range(1, 9)]
    require(workstream_ids == expected_ids, f"workstreams must be ordered exactly as {expected_ids}", errors)
    require(len(set(workstream_ids)) == 8, "workstream ids must be unique", errors)

    obligations = record["open_obligations"]
    obligation_by_id = {item["id"]: item for item in obligations}
    require(len(obligation_by_id) == len(obligations), "obligation ids must be unique", errors)

    all_gate_ids: set[str] = set()
    all_deliverable_ids: set[str] = set()
    required_gate_total = 0
    open_required_gate_total = 0

    for workstream in workstreams:
        wsid = workstream["id"]
        gates = workstream["acceptance_gates"]
        gate_ids = {gate["id"] for gate in gates}
        require(len(gate_ids) == len(gates), f"{wsid} gate ids must be unique", errors)
        require(not all_gate_ids.intersection(gate_ids), f"{wsid} reuses a gate id", errors)
        all_gate_ids.update(gate_ids)

        deliverable_ids = {item["id"] for item in workstream["deliverables"]}
        require(len(deliverable_ids) == len(workstream["deliverables"]), f"{wsid} deliverable ids must be unique", errors)
        require(not all_deliverable_ids.intersection(deliverable_ids), f"{wsid} reuses a deliverable id", errors)
        all_deliverable_ids.update(deliverable_ids)

        for deliverable in workstream["deliverables"]:
            for gate_id in deliverable["acceptance_gate_ids"]:
                require(gate_id in gate_ids, f"{deliverable['id']} references non-local gate {gate_id}", errors)

        required = sum(1 for gate in gates if gate["required"])
        classified = sum(1 for gate in gates if gate["claim_status"] in {"unknown", "supported", "refuted", "not_applicable"})
        supported = sum(1 for gate in gates if gate["claim_status"] == "supported")
        open_required = sum(1 for gate in gates if gate["required"] and gate["claim_status"] != "supported")
        required_gate_total += required
        open_required_gate_total += open_required

        expected_metrics = {
            "required_gates": required,
            "classified_gates": classified,
            "supported_gates": supported,
            "open_required_gates": open_required,
        }
        require(workstream["metrics"] == expected_metrics, f"{wsid} metrics do not match its gates", errors)

        for gate in gates:
            status = gate["claim_status"]
            evidence = gate["evidence_refs"]
            obligation_id = gate["obligation_id"]
            if status in {"supported", "refuted"}:
                require(bool(evidence), f"{gate['id']} {status} status lacks evidence", errors)
            if status == "unknown":
                require(bool(obligation_id), f"{gate['id']} unknown status lacks an obligation", errors)
            if obligation_id:
                require(obligation_id in obligation_by_id, f"{gate['id']} references missing obligation {obligation_id}", errors)
                if obligation_id in obligation_by_id:
                    require(obligation_by_id[obligation_id]["workstream_id"] == wsid, f"{gate['id']} obligation belongs to another workstream", errors)

    require(required_gate_total == 20, f"expected 20 required planning gates, got {required_gate_total}", errors)
    require(record["completion_metrics"]["obligations"]["open_required"] == open_required_gate_total, "global open-required count does not match workstream gates", errors)

    for obligation in obligations:
        require(obligation["workstream_id"] in workstream_ids, f"{obligation['id']} has unknown workstream", errors)
        require(obligation["status"] != "resolved", f"resolved obligation {obligation['id']} must not remain in open_obligations", errors)

    dependency_pairs = {(item["from"], item["to"], item["relation"]) for item in record["dependencies"]}
    require(len(dependency_pairs) == len(record["dependencies"]), "dependencies must be unique", errors)
    blocking_edges: list[tuple[str, str]] = []
    for dependency in record["dependencies"]:
        source = dependency["from"]
        target = dependency["to"]
        require(source != target, f"self-dependency {source}", errors)
        require(source in workstream_ids and target in workstream_ids, f"dependency has unknown endpoint: {source}->{target}", errors)
        if dependency["relation"] in {"requires", "blocks"}:
            blocking_edges.append((source, target))
    check_acyclic(blocking_edges, errors)

    for workstream in workstreams:
        declared = set(workstream["depends_on"])
        listed = {
            dependency["to"]
            for dependency in record["dependencies"]
            if dependency["from"] == workstream["id"] and dependency["relation"] == "requires"
        }
        require(declared == listed, f"{workstream['id']} depends_on and requires edges differ", errors)

    route_lanes = [route["lane"] for route in record["attack_routes"]]
    require(route_lanes == ["A", "B", "C", "D", "E"], "mock must contain one ordered skeleton record for lanes A-E", errors)
    require(all(route["claim_status"] == "unknown" for route in record["attack_routes"]), "mock attack routes cannot be promoted", errors)
    require(record["subjects"] == [], "planning mock cannot contain authoritative audit subjects", errors)
    require(record["assessments"] == [], "planning mock cannot contain authoritative assessments", errors)
    require(record["inventory_freeze"]["status"] == "not_frozen", "planning mock cannot freeze the authoritative denominator", errors)

    for source in record["source_artifacts"]:
        path = REPO / source["path"]
        require(path.is_file(), f"source artifact missing: {source['path']}", errors)
        if source["status"] == "pinned":
            revision = source["revision"] or ""
            require(len(revision) == 40 and all(character in "0123456789abcdef" for character in revision), f"pinned source lacks a 40-hex revision: {source['id']}", errors)
            require(bool(source["sha256"]), f"pinned source lacks a digest: {source['id']}", errors)
        if source["status"] == "current_worktree":
            require(source["revision"] is None, f"current-worktree source must not claim a committed revision: {source['id']}", errors)
        if source["sha256"] and path.is_file():
            require(sha256_file(path) == source["sha256"], f"source artifact hash mismatch: {source['path']}", errors)

    for key in ("plan_path", "pdf_path", "diagram_path", "validator_path"):
        path = REPO / record["report"][key]
        require(path.is_file(), f"report artifact missing: {record['report'][key]}", errors)

    derived = record["completion_metrics"]["derived"]
    require(not any(derived.values()), "planning mock completion flags must all remain false", errors)

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(f"ECDLP audit plan verification failed with {len(errors)} error(s).", file=sys.stderr)
        return 1

    print("ECDLP audit plan verification: PASS")
    print(f"plan_uid={record['plan_uid']}")
    print(f"workstreams={len(workstreams)} required_gates={required_gate_total} open_required={open_required_gate_total}")
    print("authoritative=false promotion_allowed=false subjects=0 assessments=0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
