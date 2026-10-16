#!/usr/bin/env python3
"""Freeze a source-only, known-answer 3-small/2-residual N83 control."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PRIOR = HERE.parent / "mixed-orbit-search-15000-20261012"
sys.path.insert(0, str(STUDY))
from chain_search_supervisor import (  # noqa: E402
    PRIMARY_ORDER, parse_point, point_add, point_mul,
)

SMALL_COUNT = 196212
ENVELOPE_COUNT = 339968
ORBIT_SIZE = 166
SMALL_INDICES = (0, 166, 332)
RESIDUAL_INDICES = (196212, 196378)
CURVE_SLUG = "icv1-f2m83-tm6151469093347-debefd74"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def collect(small_path: Path, envelope_path: Path) -> dict:
    receipts = json.loads((PRIOR / "evidence-01/objects.json").read_text())
    receipt_by_label = {row["label"]: row for row in receipts["objects"]}
    require(set(receipt_by_label) == {"small", "envelope"}, "object receipts")
    for label, path in (("small", small_path), ("envelope", envelope_path)):
        receipt = receipt_by_label[label]
        require(path.stat().st_size == receipt["size_bytes"], f"{label} byte count")
        require(sha256(path) == receipt["sha256"], f"{label} SHA-256")

    prior_audit_path = PRIOR / "audit.json"
    prior_audit = json.loads(prior_audit_path.read_text())
    require(prior_audit["status"] == "PASS"
            and prior_audit["cumulative_conservative_charge_seconds"] == 15000
            and prior_audit["small_object_sha256"] == receipt_by_label["small"]["sha256"]
            and prior_audit["envelope_object_sha256"] == receipt_by_label["envelope"]["sha256"],
            "prior search audit binding")

    selected = {}
    with gzip.open(small_path, "rt") as small, gzip.open(envelope_path, "rt") as envelope:
        small_header = json.loads(next(small))
        envelope_header = json.loads(next(envelope))
        for header, count, columns in ((small_header, SMALL_COUNT, 1182),
                                        (envelope_header, ENVELOPE_COUNT, 2048)):
            require(header["point_count"] == count
                    and header["orbit_columns"] == columns
                    and header["curve_slug"] == CURVE_SLUG
                    and header["subgroup_order"] == str(PRIMARY_ORDER)
                    and header["closure"] == "signed_frobenius", "archive header")
        for index in range(ENVELOPE_COUNT):
            envelope_line = next(envelope, None)
            require(envelope_line is not None, "short envelope archive")
            if index < SMALL_COUNT:
                small_line = next(small, None)
                require(small_line == envelope_line, "nested archive prefix differs")
            if index in SMALL_INDICES or index in RESIDUAL_INDICES:
                row = json.loads(envelope_line)
                require(row["column"] == index // ORBIT_SIZE
                        and row["phase"] == 0 and row["negative"] is False
                        and row["coefficient"] == "1", "selected orbit representative")
                point = parse_point(row["point"])
                require(point is not None and point_mul(point, PRIMARY_ORDER) is None,
                        "selected point outside primary subgroup")
                selected[index] = row
        require(next(small, None) is None and next(envelope, None) is None,
                "archive contains extra point rows")
    require(set(selected) == set(SMALL_INDICES + RESIDUAL_INDICES), "missing selected rows")
    require(len({tuple(row["point"]) for row in selected.values()}) == 5,
            "selected points overlap")

    target = None
    terms = []
    for index in SMALL_INDICES + RESIDUAL_INDICES:
        row = selected[index]
        target = point_add(target, parse_point(row["point"]))
        terms.append({"domain": "small" if index < SMALL_COUNT else "residual",
                      "global_point_index": index, **row})
    require(target is not None and point_mul(target, PRIMARY_ORDER) is None,
            "target identity or outside primary subgroup")
    target_words = [str(word) for word in target]
    require(parse_point(target_words) == target, "target outside curve")

    return {
        "schema": "n83.mixed-orbit-known-witness/v1",
        "scope": "known_answer_control_only",
        "curve_slug": CURVE_SLUG,
        "subgroup_order": str(PRIMARY_ORDER),
        "archive": {
            label: {"s3_uri": receipt_by_label[label]["s3_uri"],
                    "sha256": receipt_by_label[label]["sha256"],
                    "size_bytes": receipt_by_label[label]["size_bytes"]}
            for label in ("small", "envelope")
        },
        "prior_search_audit_sha256": sha256(prior_audit_path),
        "equation": {"a": "0", "b": "1", "target": target_words,
                     "small_summands": 3, "residual_summands": 2},
        "terms": terms,
        "verification": {
            "nested_point_prefix_equal": True,
            "selected_columns_distinct": True,
            "selected_points_on_curve_and_in_subgroup": True,
            "nonidentity_target_on_curve_and_in_subgroup": True,
            "ordered_selected_point_sum_equals_target": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, required=True)
    parser.add_argument("--envelope", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, default=HERE / "fixture.json")
    parser.add_argument("--check", action="store_true", help="compare without writing")
    args = parser.parse_args()
    rendered = json.dumps(collect(args.small, args.envelope), indent=2,
                          sort_keys=True) + "\n"
    if args.check:
        require(args.fixture.read_text() == rendered, "frozen fixture differs")
        print("PASS: archive-bound known-answer fixture")
    else:
        args.fixture.write_text(rendered)
        print(f"Wrote {args.fixture}")


if __name__ == "__main__":
    main()
