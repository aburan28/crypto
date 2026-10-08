"""Canonical JSON, hashing and the identifiers ICMS mints.

Canonical JSON follows the convention both repositories already use for
identity hashing: sorted keys, compact separators, ASCII, no floats in hashed
inputs (a float in an identity is refused, because its decimal rendering is
not unique across languages).
"""
from __future__ import annotations

import hashlib
import json
from typing import Any


class CanonicalError(ValueError):
    pass


def _check(obj: Any, path: str = "$") -> None:
    if isinstance(obj, float):
        raise CanonicalError(f"{path}: floats are not allowed in hashed identity inputs; use an exact decimal string")
    if isinstance(obj, dict):
        for k, v in obj.items():
            if not isinstance(k, str):
                raise CanonicalError(f"{path}: non-string key {k!r}")
            _check(v, f"{path}.{k}")
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            _check(v, f"{path}[{i}]")


def canonical_json(obj: Any, strict: bool = True) -> str:
    if strict:
        _check(obj)
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_hex(obj: Any, strict: bool = True) -> str:
    return hashlib.sha256(canonical_json(obj, strict).encode("ascii")).hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str | None:
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return None
    return h.hexdigest()


# --- identifiers -------------------------------------------------------------
#
# ICS1h<12>   spec id: the canonical spec minus its non-identity fields
# ENV1h<12>   environment class: the stable host facts a wall-clock pair shares
# W<12>       workload id: curve + targets + law + seeds + cold target count,
#             formed exactly as cryptanalysis AGENTS.md prescribes
# IC1...h<12> candidate id: cryptanalysis AGENTS.md's compact label, hashed
#             over the canonical candidate record (completed after the run,
#             because fb<B> is the measured usable point count)
# ICR1h<12>   run record id (content hash of the record minus its own id)

def short(hexdigest: str, n: int = 12) -> str:
    return hexdigest[:n]


def spec_id(identity_view: dict) -> str:
    return "ICS1h" + short(sha256_hex(identity_view))


def env_class_id(env_class: dict) -> str:
    return "ENV1h" + short(sha256_hex(env_class, strict=False))


def workload_id(workload_record: dict) -> str:
    return "W" + short(sha256_hex(workload_record))


def record_id(record_without_id: dict) -> str:
    return "ICR1h" + short(sha256_hex(record_without_id, strict=False))


def candidate_label(n: int | None, curve_tag: str, fb_points: int | None, m: int, pdp_code: str,
                    rc_code: str, la_code: str, td_code: str, iso: int, candidate_record: dict) -> str:
    """IC1N<n>C<tag>fb<B>PDP<m><solver>RC<rc>LA<la>TD<td>ISO<0|1>h<12hex>.

    ``n`` is the field degree (log2 of the field size), never a subgroup size.
    ``fb_points`` is the actual usable point count; ``None`` gives no label,
    because a nominal dimension is never an fb count (cryptanalysis AGENTS.md).
    """
    if n is None or fb_points is None:
        return None
    digest = short(sha256_hex(candidate_record))
    return f"IC1N{n}C{curve_tag}fb{fb_points}PDP{m}{pdp_code}RC{rc_code}LA{la_code}TD{td_code}ISO{iso}h{digest}"
