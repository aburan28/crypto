#!/usr/bin/env python3
"""Compare distinguished-point corpora without misframing versioned files.

Benchmark jobs write ``dp-<variant>.bin`` under one result directory.  Legacy
sigma corpora are headerless 32-byte records.  Witness-v2 and table-walk-v3
corpora carry a 16-byte header and use the record size declared there.  The
comparison is over sorted complete records after that header.
"""

from __future__ import annotations

import hashlib
import pathlib
import struct
import sys


HEADER = struct.Struct("<8sII")
FORMATS = {
    b"ECC2KDP2": (2, 72, "v2"),
    b"ECC2KDT3": (3, 32, "table-v3"),
}


def corpus_records(path: pathlib.Path) -> tuple[list[bytes], str, int]:
    """Return complete records, format name and payload offset for *path*."""
    data = path.read_bytes()
    offset, stride, format_name = 0, 32, "v1"
    if len(data) >= 8 and data[:8] in FORMATS:
        if len(data) < HEADER.size:
            raise ValueError(f"{path}: truncated corpus header")
        magic, version, declared_stride = HEADER.unpack_from(data)
        expected_version, expected_stride, format_name = FORMATS[magic]
        if version != expected_version or declared_stride != expected_stride:
            raise ValueError(
                f"{path}: invalid {format_name} header: version {version}, "
                f"record bytes {declared_stride}"
            )
        offset, stride = HEADER.size, declared_stride
    payload = data[offset:]
    if len(payload) % stride:
        raise ValueError(
            f"{path}: {len(payload)} payload bytes are not a multiple of {stride}"
        )
    return [payload[i : i + stride] for i in range(0, len(payload), stride)], format_name, offset


def corpus_identity(path: pathlib.Path) -> dict[str, object]:
    records, format_name, offset = corpus_records(path)
    digest = hashlib.sha256(b"".join(sorted(records))).hexdigest()
    return {
        "records": len(records),
        "sha256": digest,
        "format": format_name,
        "headerBytes": offset,
    }


def compare(result_dir: pathlib.Path, variants: list[str]) -> tuple[list[str], bool]:
    if not variants:
        raise ValueError("at least one variant is required")
    rows: list[str] = []
    reference = None
    valid = True
    for name in variants:
        path = result_dir / f"dp-{name}.bin"
        if not path.is_file():
            rows.append(f"{name:<20} missing")
            valid = False
            continue
        identity = corpus_identity(path)
        if reference is None:
            reference = identity
        identical = identity == reference
        valid = valid and identical
        rows.append(
            "%-20s %8d records  sha256 %s  %s"
            % (
                name,
                identity["records"],
                str(identity["sha256"])[:16],
                "IDENTICAL to %s" % variants[0]
                if identical
                else "DIFFERS from %s" % variants[0],
            )
        )
    return rows, valid


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print("usage: dp_identity.py RESULT_DIR VARIANT [VARIANT ...]", file=sys.stderr)
        return 2
    try:
        rows, valid = compare(pathlib.Path(argv[1]), argv[2:])
    except (OSError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print("\n".join(rows))
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
