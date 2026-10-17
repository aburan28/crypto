#!/usr/bin/env python3
"""Extract exact full-orbit domains for the ordered N83 encoding control."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
FIXTURE = STUDY / "verification/mixed-orbit-known-witness-source-20261012/fixture.json"
sys.path.insert(0, str(STUDY))
from chain_search_supervisor import parse_point, point_add  # noqa: E402

ORBIT_SIZE = 166
SMALL_ORBITS = 1182
ENVELOPE_ORBITS = 2048


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect(small_path, envelope_path):
    fixture_bytes = FIXTURE.read_bytes()
    fixture = json.loads(fixture_bytes)
    for label, path in (("small", small_path), ("envelope", envelope_path)):
        expected = fixture["archive"][label]
        require(path.stat().st_size == expected["size_bytes"], f"{label} compressed length")
        require(sha256(path) == expected["sha256"], f"{label} compressed SHA-256")
    small_reps = []
    residual_reps = []
    selected = {}
    with gzip.open(small_path, "rt") as small, gzip.open(envelope_path, "rt") as envelope:
        small_header = json.loads(next(small))
        envelope_header = json.loads(next(envelope))
        for header, orbits in ((small_header, SMALL_ORBITS),
                               (envelope_header, ENVELOPE_ORBITS)):
            require(header["point_count"] == orbits * ORBIT_SIZE
                    and header["orbit_columns"] == orbits
                    and header["curve_slug"] == fixture["curve_slug"]
                    and header["closure"] == "signed_frobenius", "archive header")
        for index in range(ENVELOPE_ORBITS * ORBIT_SIZE):
            row_bytes = next(envelope, None)
            require(row_bytes is not None, "short envelope archive")
            if index < SMALL_ORBITS * ORBIT_SIZE:
                require(next(small, None) == row_bytes, "non-nested small archive")
            if index % ORBIT_SIZE == 0:
                row = json.loads(row_bytes)
                require(row["column"] == index // ORBIT_SIZE
                        and row["phase"] == 0 and row["negative"] is False,
                        "orbit representative label")
                x = row["point"][0]
                if index < SMALL_ORBITS * ORBIT_SIZE:
                    small_reps.append(x)
                else:
                    residual_reps.append(x)
                if index in (0, 166, 332, 196212, 196378):
                    selected[index] = row
        require(next(small, None) is None and next(envelope, None) is None,
                "extra archive point rows")
    require(len(small_reps) == SMALL_ORBITS
            and len(residual_reps) == ENVELOPE_ORBITS - SMALL_ORBITS,
            "representative count")
    require(len(set(small_reps)) == SMALL_ORBITS
            and len(set(residual_reps)) == ENVELOPE_ORBITS - SMALL_ORBITS
            and set(small_reps).isdisjoint(residual_reps), "representatives overlap")
    for term in fixture["terms"]:
        index = term["global_point_index"]
        require(selected[index] == {key: term[key] for key in selected[index]},
                f"selected archive term {index}")
    ordered_terms = [fixture["terms"][i] for i in (1, 0, 2, 3, 4)]
    require(int(ordered_terms[0]["point"][0]) < int(ordered_terms[1]["point"][0]),
            "first-pair order")
    partial = None
    partial_xs = []
    for term in ordered_terms:
        partial = point_add(partial, parse_point(term["point"]))
        require(partial is not None, "intermediate identity")
        partial_xs.append(str(partial[0]))
    require([str(v) for v in partial] == fixture["equation"]["target"], "target sum")
    return {
        "schema": "n83.fixed-witness-domains/v1",
        "scope": "known_answer_encoding_control_only",
        "identity_mask": 0,
        "fixture_sha256": hashlib.sha256(fixture_bytes).hexdigest(),
        "archive_sha256": {"small": sha256(small_path),
                           "envelope": sha256(envelope_path)},
        "ordered_indices": [term["global_point_index"] for term in ordered_terms],
        "ordered_xs": [term["point"][0] for term in ordered_terms],
        "partial_xs": partial_xs[1:-1],
        "target_x": fixture["equation"]["target"][0],
        "small_reps": small_reps,
        "residual_reps": residual_reps,
        "natural_relations": None,
        "selected_best_total_runtime": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, required=True)
    parser.add_argument("--envelope", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "domains.json")
    args = parser.parse_args()
    output = collect(args.small, args.envelope)
    with args.output.open("x") as stream:
        json.dump(output, stream, sort_keys=True, indent=2)
        stream.write("\n")
    print(json.dumps({"status": "PASS", "small_orbits": len(output["small_reps"]),
                      "residual_orbits": len(output["residual_reps"]),
                      "output_sha256": sha256(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
