#!/usr/bin/env python3
"""Freeze #747 bases and a fresh, disjoint point-only target stream."""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
SWEEP = PARENT / "compact_base_sweep_20260925"
DOMAIN = b"ECC2K-COMPACT-FOUR-SUM-ORACLE-20260925-v1/"
ARMS = {
    (37, 3): ("train-n37-R3-c512", 222,
              "2722f3c7271ea410e47ff9d20b67496a9a28120804f8b790e0729d8c86ea6ddb"),
    (41, 8): ("train-n41-R8-c4096", 656,
              "cbdd871504ebac2f36e16ce652bc6dfe9623d3fd62723b05f7112d3d94614520"),
    (41, 12): ("train-n41-R12-c512", 984,
               "a3856771a5ea1ad3a1e2e262eebcc1b2ccde8293e53af4360268cdcf72ab7fa8"),
}
ORDERS = {37: 230603167, 41: 549756390943}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def once(path: Path, data: bytes) -> None:
    if path.exists():
        assert path.read_bytes() == data, f"frozen file changed: {path}"
    else:
        path.write_bytes(data)


def compact(value) -> bytes:
    return json.dumps(value, separators=(",", ":")).encode()


def main() -> None:
    # Use exactly the independently checked affine arithmetic from #737.
    # An import would create pycache in the archived source folder, so exec it.
    source = SWEEP / "evidence/source/pr737_independent_math.py"
    namespace = {"__name__": "pr737_independent_math", "__file__": str(source)}
    exec(compile(source.read_bytes(), str(source), "exec"), namespace)
    curve_class = namespace["Curve"]

    headers = {}
    for (n, r), (run_name, size, expected_hash) in ARMS.items():
        archive = SWEEP / "evidence/runs" / run_name / "producer.stdout.jsonl.gz"
        with gzip.open(archive, "rt") as stream:
            result = json.loads(stream.readline())
        header = result["compact_orbit_base_header"]
        assert (header["n"], header["a"], header["orbit_columns"],
                header["factor_base_points"]) == (n, 0, r, size)
        assert int(header["subgroup_order"]) == ORDERS[n]
        assert sha(compact(header["factor_base_point_coordinates"])) == expected_hash
        assert len(set(map(tuple, header["factor_base_point_coordinates"]))) == size
        once(HERE / f"base_n{n}_R{r}.json", compact(header) + b"\n")
        headers[(n, r)] = header

    spec = json.loads((SWEEP / "input_spec.json").read_bytes())
    manifest = {"schema_version": "1.0", "domain": DOMAIN.decode(),
                "source_pr747_commit": "fc27150df3238b6863ed5618c721e7fd8b6ce403",
                "input_spec_sha256": sha((SWEEP / "input_spec.json").read_bytes()),
                "independent_math_sha256": sha(source.read_bytes()),
                "arms": {}}
    for n, order in ORDERS.items():
        header = headers[(n, 3 if n == 37 else 8)]
        curve = curve_class(n, 0, header["field_modulus_low_terms"])
        generator = tuple(header["generator"])
        assert curve.is_on_curve(generator)
        assert curve.scalar(generator, order) is None
        exclusions = set(map(int, (SWEEP / f"training_scalars_n{n}.txt").read_text().split()))
        exclusions.update(int(row["scalar_validator_only"])
                          for row in spec["holdouts"][str(n)])
        assert len(exclusions) == 4099
        seen = set(exclusions)
        scalars = []
        counter = 0
        while len(scalars) < 512:
            digest = hashlib.sha256(DOMAIN + str(n).encode() + b"/"
                                    + str(counter).encode()).digest()
            scalar = 1 + int.from_bytes(digest, "big") % (order - 1)
            counter += 1
            if scalar in seen:
                continue
            seen.add(scalar)
            scalars.append(scalar)
        points = [curve.scalar(generator, scalar) for scalar in scalars]
        assert all(point is not None and curve.is_on_curve(point) for point in points)
        assert len(set(points)) == 512
        scalar_data = "".join(f"{scalar}\n" for scalar in scalars).encode()
        point_data = "".join(compact(point).decode() + "\n" for point in points).encode()
        once(HERE / f"target_scalars_n{n}.txt", scalar_data)
        once(HERE / f"target_points_n{n}.jsonl", point_data)
        manifest["arms"][str(n)] = {
            "generator": generator, "subgroup_order": order,
            "hash_counter_exclusive": counter,
            "target_count": 512,
            "target_scalars_sha256": sha(scalar_data),
            "target_points_sha256": sha(point_data),
        }
    # The R12 extractor consumes its whole point file; freeze precisely the
    # preregistered first 128 lines as a separate, byte-identical prefix.
    n41_lines = (HERE / "target_points_n41.jsonl").read_bytes().splitlines(keepends=True)
    assert len(n41_lines) == 512
    once(HERE / "target_points_n41_R12.jsonl", b"".join(n41_lines[:128]))
    manifest_data = json.dumps(manifest, sort_keys=True, indent=2).encode() + b"\n"
    once(HERE / "input_manifest.json", manifest_data)
    print(json.dumps({"input_manifest_sha256": sha(manifest_data),
                      "arms": manifest["arms"]}, sort_keys=True))


if __name__ == "__main__":
    main()
