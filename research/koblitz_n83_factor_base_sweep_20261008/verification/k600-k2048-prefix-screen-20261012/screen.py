#!/usr/bin/env python3
"""Verify the retained K600/K2048 prefix and count its mixed support."""

import argparse
from decimal import Decimal, localcontext
import gzip
import hashlib
import json
from math import comb
from pathlib import Path


HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
PILOT_MANIFEST = STUDY / "pilot-01/manifest.json"
PILOT_UPLOAD = STUDY / "pilot-01/upload-receipt.json"
ENVELOPE_MANIFEST = STUDY / "verification/remaining-budget-k2048-20261010/manifest.json"
PRIOR_AUDIT = STUDY / "verification/mixed-orbit-search-15000-20261012/audit.json"
S3_HEAD = HERE / "s3-head.json"
SMALL_SHA256 = "3e85c9e419fcddd3c9faed3c84470bf96fcad9178e48e1d73629369335841ecf"
ENVELOPE_SHA256 = "aacb62784dee0c9f5e06328dbb85a4b6a1c7243e778328da82900e7284f3a74a"
SMALL_BLAKE3 = "b4dda45e6c3edfca673db6f3f82c16cce8c860f5e62c092155258df399a0363e"
ENVELOPE_BLAKE3 = "82990cda5d57b1c4c4b635673a881ad233c24f707e7ac644a11478aeddddf851"
PRIOR_AUDIT_SHA256 = "32420a62bb6ccf8f9df7f03dea9040e518606740968110dc8c0ac47cb50cafd5"
SMALL_COLUMNS, ENVELOPE_COLUMNS, ORBIT_SIZE = 600, 2048, 166
ORDER = 2417851639230796216685689


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sole_base(manifest: dict, columns: int) -> dict:
    rows = [row for row in manifest["bases"]
            if (row["a"], row["policy"], row["seed"], row["columns"])
            == (0, "public_x_hash", 2026100801, columns)]
    require(len(rows) == 1, f"unique K={columns} manifest row")
    return rows[0]


def count_support(small_points: int, residual_points: int) -> list[dict]:
    rows = []
    for summands in (5, 6):
        for residuals in (0, 1, 2):
            count = (comb(small_points + summands - residuals - 1,
                          summands - residuals)
                     * comb(residual_points + residuals - 1, residuals))
            with localcontext() as context:
                context.prec = 30
                ratio = f"{Decimal(count) / Decimal(ORDER - 1):.9E}"
            rows.append({
                "summands": summands, "small_summands": summands - residuals,
                "residual_summands": residuals,
                "unordered_point_multisets_with_repetition": str(count),
                "count_divided_by_nonidentity_subgroup_points": ratio,
                "uniform_nonidentity_hit_probability_ceiling": {
                    "numerator": str(min(count, ORDER - 1)),
                    "denominator": str(ORDER - 1),
                },
            })
    return rows


def screen(small: Path, envelope: Path) -> dict:
    pilot = json.loads(PILOT_MANIFEST.read_text())
    upload = json.loads(PILOT_UPLOAD.read_text())
    large = json.loads(ENVELOPE_MANIFEST.read_text())
    audit = json.loads(PRIOR_AUDIT.read_text())
    head = json.loads(S3_HEAD.read_text())
    require(sha256(PRIOR_AUDIT) == PRIOR_AUDIT_SHA256
            and audit["status"] == "PASS"
            and audit["cumulative_conservative_charge_seconds"] == 15000,
            "prior N83 audit identity")
    a, b = sole_base(pilot, SMALL_COLUMNS), sole_base(large, ENVELOPE_COLUMNS)
    require(a["points"] == SMALL_COLUMNS * ORBIT_SIZE
            and a["bytes"] == small.stat().st_size == 3025798
            and a["compressed_blake3"] == SMALL_BLAKE3
            and b["points"] == ENVELOPE_COLUMNS * ORBIT_SIZE
            and b["bytes"] == envelope.stat().st_size == 10313647
            and b["compressed_blake3"] == ENVELOPE_BLAKE3,
            "manifest sizes or compressed-object names")
    require(sha256(small) == SMALL_SHA256 and sha256(envelope) == ENVELOPE_SHA256,
            "downloaded compressed bytes")
    require(head == {"ContentLength": 3025798,
                     "ETag": '"359c027f8c0be18a8e41696f6c5c7aed"'},
            "observed S3 head differs")
    receipts = [row for row in upload["objects"] if row["s3_uri"] == a["s3_uri"]]
    require(upload["status"] == "PASS" and len(receipts) == 1
            and receipts[0]["compressed_blake3"] == SMALL_BLAKE3
            and receipts[0]["downloaded_hash_matches"] is True,
            "original K600 S3 round-trip receipt")

    with gzip.open(small, "rb") as left, gzip.open(envelope, "rb") as right:
        h_small, h_envelope = json.loads(next(left)), json.loads(next(right))
        for key, value in (
            ("n", 83), ("curve_a", 0), ("curve_slug", "icv1-f2m83-tm6151469093347-debefd74"),
            ("policy", "public_x_hash"), ("seed", 2026100801),
            ("closure", "signed_frobenius"), ("subgroup_order", str(ORDER)),
        ):
            require(h_small[key] == h_envelope[key] == value, f"header {key}")
        require(h_small["point_count"] == a["points"]
                and h_small["orbit_columns"] == SMALL_COLUMNS
                and h_envelope["point_count"] == b["points"]
                and h_envelope["orbit_columns"] == ENVELOPE_COLUMNS,
                "header counts")
        for key in ("accepted_candidate_indices", "representatives",
                    "source_points_before_cofactor_clearing"):
            require(h_small[key] == h_envelope[key][:len(h_small[key])],
                    f"ordered header {key} prefix")
        for index in range(a["points"]):
            require(next(left, None) == next(right, None), f"point row {index} differs")
        require(next(left, None) is None, "K600 archive has extra rows")
        for index in range(a["points"], b["points"]):
            require(next(right, None) is not None, f"K2048 archive short at {index}")
        require(next(right, None) is None, "K2048 archive has extra rows")

    residual_points = b["points"] - a["points"]
    return {
        "schema": "n83.k600-k2048-exact-prefix-screen/v1",
        "scope": "source_only_nested_pair_and_support",
        "small_s3_uri": a["s3_uri"],
        "envelope_s3_uri": b["s3_uri"],
        "small_object_sha256": SMALL_SHA256,
        "envelope_object_sha256": ENVELOPE_SHA256,
        "pilot_manifest_sha256": sha256(PILOT_MANIFEST),
        "pilot_upload_receipt_sha256": sha256(PILOT_UPLOAD),
        "envelope_manifest_sha256": sha256(ENVELOPE_MANIFEST),
        "prior_search_audit_sha256": PRIOR_AUDIT_SHA256,
        "s3_head_receipt_sha256": sha256(S3_HEAD),
        "small_columns": SMALL_COLUMNS,
        "small_points": a["points"],
        "residual_columns": ENVELOPE_COLUMNS - SMALL_COLUMNS,
        "residual_points": residual_points,
        "ordered_header_and_all_point_rows_prefix": True,
        "support": count_support(a["points"], residual_points),
        "solver_search_executed": False,
        "selected_best_total_runtime": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--small", type=Path, required=True)
    parser.add_argument("--envelope", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "screen.json")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(screen(args.small, args.envelope), indent=2,
                          sort_keys=True) + "\n"
    if args.check:
        require(args.output.read_text() == rendered, "saved screen differs")
        print("PASS: exact K600/K2048 prefix and support")
    else:
        args.output.write_text(rendered)
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
