#!/usr/bin/env python3
"""Bind both primary SEC 2 editions to the existing exact curve/trait catalog.

Requires pdftotext (Poppler). This script performs parameter comparisons and
copies existing trait records; it performs no curve search or factoring.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = {
    "2000_v1": (HERE / "sources/sec2-v1.pdf", "https://www.secg.org/SEC2-Ver-1.0.pdf"),
    "2010_v2": (HERE / "sources/sec2-v2.pdf", "https://www.secg.org/sec2-v2.pdf"),
}
SECTION = re.compile(
    r"^\s*([23]\.\d+\.\d+)\s+Recommended Parameters\s+"
    r"(sec[pt]\d{3}[kr]\d)\s*$",
    re.MULTILINE,
)
HEX_LINE = re.compile(r"^\s*[0-9A-Fa-f ]+\s*$")
TRIAL_LIMIT = 100_000


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_model(curve: dict) -> str:
    return hashlib.sha256(curve["model_json"].encode("utf-8")).hexdigest()


def as_int(value: object) -> int:
    if isinstance(value, int):
        return value
    return int(str(value), 0) if str(value).startswith("0x") else int(str(value))


def field_values(section: str, label: str) -> list[str]:
    lines = section.splitlines()
    matches: list[str] = []
    start = re.compile(r"^\s*" + re.escape(label) + r"\s*=\s*([0-9A-Fa-f ]+)\s*$")
    for i, line in enumerate(lines):
        m = start.match(line)
        if not m:
            continue
        chunks = [m.group(1).replace(" ", "")]
        for tail in lines[i + 1 :]:
            if not HEX_LINE.fullmatch(tail) or not tail.strip():
                break
            chunks.append(tail.strip().replace(" ", ""))
        matches.append("".join(chunks).lower())
    return matches


def reduction_exponents(section: str) -> list[int] | None:
    m = re.search(r"^\s*f\s*\(x\)\s*=\s*(.+)$", section, re.MULTILINE)
    if not m:
        return None
    expression = m.group(1)
    exponents = [int(x) if x else 1 for x in re.findall(r"x\s*(\d*)", expression)]
    if re.search(r"(?:^|\+)\s*1\s*$", expression):
        exponents.append(0)
    return sorted(exponents, reverse=True)


def pdf_sections(path: Path) -> dict[str, dict]:
    text = subprocess.check_output(["pdftotext", "-layout", str(path), "-"], text=True)
    headings = list(SECTION.finditer(text))
    result: dict[str, dict] = {}
    for i, heading in enumerate(headings):
        name = heading.group(2)
        if name in result:
            raise ValueError(f"{path}: duplicate parameter section {name}")
        body = text[heading.end() : headings[i + 1].start() if i + 1 < len(headings) else len(text)]
        result[name] = {
            "section": heading.group(1),
            "p": field_values(body, "p"),
            "a": field_values(body, "a"),
            "b": field_values(body, "b"),
            "G": field_values(body, "G"),
            "n": field_values(body, "n"),
            "h": field_values(body, "h"),
            "seed": field_values(body, "S"),
            "reduction_exponents": reduction_exponents(body),
        }
    return result


def only(record: dict, key: str) -> str:
    values = record[key]
    if len(values) != 1:
        raise ValueError(f"{key}: expected one value, got {values}")
    return values[0]


def verify_record(name: str, source: dict, registry: dict) -> dict:
    representation = registry["representations"][0]
    curve = representation["curve"]
    field = representation["field"]
    bits = (as_int(field["degree"]) if name.startswith("sect")
            else as_int(field["characteristic"]).bit_length())
    octets = (bits + 7) // 8
    encoded = source["G"][-1]
    if not encoded.startswith("04") or len(encoded) != 2 + 4 * octets:
        raise ValueError(f"{name}: invalid SEC 2 uncompressed generator length {len(encoded)}")
    x = int(encoded[2 : 2 + 2 * octets], 16)
    y = int(encoded[2 + 2 * octets :], 16)
    coefficient_a = curve["coefficients"][1 if name.startswith("sect") else 3]
    coefficient_b = curve["coefficients"][4]
    checks = {
        "coefficient_a": int(only(source, "a"), 16) == as_int(coefficient_a),
        "coefficient_b": int(only(source, "b"), 16) == as_int(coefficient_b),
        "generator_x": x == int(curve["generator"][0], 16),
        "generator_y": y == int(curve["generator"][1], 16),
        "subgroup_order": int(only(source, "n"), 16) == as_int(curve["subgroup_order"]),
        "cofactor": int(only(source, "h"), 16) == as_int(curve["cofactor"]),
    }
    if name.startswith("secp"):
        checks["field_prime"] = int(only(source, "p"), 16) == as_int(registry["params"]["p"])
    else:
        checks["field_degree"] = bits == as_int(field["degree"])
        checks["reduction_polynomial"] = source["reduction_exponents"] == field["modulus_exponents"]
    if not all(checks.values()):
        raise ValueError(f"{name}: primary-source mismatch {checks}")
    return checks


def primes_through(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (((limit - p * p) // p) + 1)
    return [p for p, is_prime in enumerate(sieve) if is_prime]


def cheon_advisory(n: int, trial_primes: list[int]) -> dict:
    """A bounded check of SEC 1 v2's r > n^(19/20) recommendation.

    If the residual after exact trial division is below the threshold, no
    unexamined prime factor can pass it. A larger residual is unresolved.
    """
    result = {}
    for label, value in (("n_minus_1", n - 1), ("n_plus_1", n + 1)):
        residual = value
        factors = []
        for p in trial_primes:
            if residual == 1:
                break
            power = 0
            while residual % p == 0:
                residual //= p
                power += 1
            if power:
                factors.append([str(p), power])
        fails = residual**20 <= n**19
        result[label] = {
            "status": "fails_advisory_bound" if fails else "undetermined",
            "trial_division_limit": TRIAL_LIMIT,
            "small_prime_factors": factors,
            "residual": str(residual),
            "residual_bits": residual.bit_length(),
            "logic": "Every remaining prime factor is at most the residual; no assertion that the residual is prime.",
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="check the committed inventory without writing")
    args = parser.parse_args()
    registry_path = ROOT / "docs/curves/registry.json"
    traits_path = ROOT / "docs/curves/traits.json"
    covers_path = ROOT / "docs/curves/covers.json"
    transfers_path = ROOT / "docs/curves/jacobian-transfers.json"
    registry = json.loads(registry_path.read_text())
    traits = json.loads(traits_path.read_text())
    covers = json.loads(covers_path.read_text())
    transfers = json.loads(transfers_path.read_text())
    by_name: dict[str, list[dict]] = {}
    for curve in registry["curves"]:
        for name in curve.get("standard_names", []):
            by_name.setdefault(name, []).append(curve)
    by_slug = {curve["slug"]: curve for curve in traits["curves"]}
    cover_by_slug = {curve["slug"]: curve for curve in covers["curves"]}
    transfer_by_slug = {curve["slug"]: curve for curve in transfers["curves"]}
    primary = {edition: pdf_sections(path) for edition, (path, _) in SOURCES.items()}
    if len(primary["2000_v1"]) != 33 or len(primary["2010_v2"]) != 20:
        raise ValueError("SEC 2 primary-source section count changed")
    if not set(primary["2010_v2"]).issubset(primary["2000_v1"]):
        raise ValueError("2010 SEC 2 names are not a subset of 2000 SEC 2 names")
    trial_primes = primes_through(TRIAL_LIMIT)
    rows = []
    for name in sorted(primary["2000_v1"]):
        matches = by_name.get(name, [])
        if len(matches) != 1:
            raise ValueError(f"{name}: expected one registered curve, got {len(matches)}")
        curve = matches[0]
        trait = by_slug.get(curve["slug"])
        if trait is None:
            raise ValueError(f"{name}: no trait record for {curve['slug']}")
        if as_int(trait["order"]) != as_int(curve["order"]) or as_int(trait["trace"]) != as_int(curve["trace"]):
            raise ValueError(f"{name}: trait row differs from current curve order or trace")
        cover = cover_by_slug.get(curve["slug"])
        transfer = transfer_by_slug.get(curve["slug"])
        if cover is None or cover["model_sha256"] != sha256_model(curve):
            raise ValueError(f"{name}: missing or mismatched cover record")
        if transfer is None or transfer["icv1"] != curve["icv1"]:
            raise ValueError(f"{name}: missing or mismatched transfer record")
        source = {}
        for edition, sections in primary.items():
            if name in sections:
                record = sections[name]
                if (name[-2] == "r") != (len(record["seed"]) == 1):
                    raise ValueError(f"{name} in {edition}: unexpected seed count")
                source[edition] = {
                    "section": record["section"],
                    "seed_hex": only(record, "seed") if record["seed"] else None,
                    "checks": verify_record(name, record, curve),
                }
        if len(source) == 2 and source["2000_v1"]["seed_hex"] != source["2010_v2"]["seed_hex"]:
            raise ValueError(f"{name}: seed changed between editions")
        rows.append({
            "name": name,
            "listed_in": sorted(source),
            "source": source,
            "identity": {"icv1": curve["icv1"], "slug": curve["slug"],
                         "representations": curve["representations"]},
            "catalog": curve,
            "traits": trait,
            "cover": cover,
            "jacobian_transfer": transfer,
            "cheon_n_pm_1_advisory": cheon_advisory(as_int(curve["representations"][0]["curve"]["subgroup_order"]), trial_primes),
        })
    result = {
        "schema": "sec2-curve-inventory/1",
        "sources": {edition: {"title": "SEC 2: Recommended Elliptic Curve Domain Parameters",
                              "url": url, "sha256": sha256(path), "path": str(path.relative_to(HERE))}
                    for edition, (path, url) in SOURCES.items()},
        "catalog": {"registry_path": str(registry_path.relative_to(ROOT)),
                    "registry_sha256": sha256(registry_path),
                    "traits_path": str(traits_path.relative_to(ROOT)),
                    "traits_sha256": sha256(traits_path),
                    "traits_registry_sha256": traits["registry_sha256"],
                    "covers_path": str(covers_path.relative_to(ROOT)),
                    "covers_sha256": sha256(covers_path),
                    "covers_registry_sha256": covers["registry_sha256"],
                    "transfers_path": str(transfers_path.relative_to(ROOT)),
                    "transfers_sha256": sha256(transfers_path),
                    "transfers_registry_sha256": transfers["registry_sha256"]},
        "counts": {"2000_v1": 33, "2010_v2": 20, "historical_only": 13},
        "curves": rows,
    }
    output = HERE / "inventory.json"
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text() != serialized:
            raise ValueError(f"{output}: inventory is missing or stale")
        print(f"{output}: current")
    else:
        output.write_text(serialized)
        print(f"{output}: {len(rows)} exact curve/trait records; core model fields match both editions")


if __name__ == "__main__":
    main()
