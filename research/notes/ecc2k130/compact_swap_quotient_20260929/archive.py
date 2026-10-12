#!/usr/bin/env python3
"""Seal hosted swap-grid raw files, receipts and replay into deterministic tarballs."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import tarfile


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equivalent(left, right) -> bool:
    """Compare replay receipts exactly except for libm rounding of derived floats."""
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(
            equivalent(left[key], right[key]) for key in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(
            equivalent(a, b) for a, b in zip(left, right))
    if isinstance(left, float):
        return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12)
    return left == right


def pack(source: Path, target: Path) -> None:
    assert source.is_dir() and not target.exists()
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0,
                           compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as archive:
                for path in sorted(source.rglob("*")):
                    if not path.is_file():
                        continue
                    payload = path.read_bytes()
                    info = tarfile.TarInfo(str(Path(source.name) / path.relative_to(source)))
                    info.size = len(payload)
                    info.mtime = 0
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mode = 0o644
                    archive.addfile(info, io.BytesIO(payload))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--head", required=True)
    parser.add_argument("--checkout-merge", required=True)
    parser.add_argument("--main-parent", required=True)
    parser.add_argument("--run-url", required=True)
    parser.add_argument("--second-host", type=Path)
    parser.add_argument("--frozen", type=Path,
                        default=Path(__file__).with_name("FROZEN.json"))
    args = parser.parse_args()
    assert not args.out.exists(), "never overwrite an evidence archive"
    args.out.mkdir(parents=True)
    for folder in ("raw", "receipts", "hosts"):
        (args.out / folder).mkdir()
    if args.second_host is not None:
        (args.out / "second_host_replay").mkdir()
    cases = {}
    for n in (41, 53):
        case = f"n{n}_L1024"
        source = args.runs / case
        assert source.is_dir(), case
        receipt_src = source / "receipt.json"
        receipt = json.loads(receipt_src.read_text())
        assert receipt["status"] in ("PASS", "FAIL"), case
        if receipt["status"] == "PASS":
            assert receipt["n"] == n and receipt["L"] == 1024, case
        host_src = source / "host.json"
        host = json.loads(host_src.read_text()) if host_src.exists() else None
        if host is not None:
            assert host["git_head"] == args.checkout_merge, case
        raw = args.out / "raw" / f"{case}.tar.gz"
        pack(source, raw)
        receipt_dst = args.out / "receipts" / f"{case}.json"
        host_dst = args.out / "hosts" / f"{case}.json"
        shutil.copyfile(receipt_src, receipt_dst)
        if host is not None:
            shutil.copyfile(host_src, host_dst)
        cases[case] = {
            "status": receipt["status"],
            "raw_path": str(raw.relative_to(args.out)),
            "raw_bytes": raw.stat().st_size,
            "raw_sha256": sha(raw),
            "receipt_path": str(receipt_dst.relative_to(args.out)),
            "receipt_sha256": sha(receipt_dst),
            "host_path": str(host_dst.relative_to(args.out)) if host is not None else None,
            "host_sha256": sha(host_dst) if host is not None else None,
            "k_grid": receipt.get("k_grid"),
            "targets_verified_per_arm": receipt.get("targets_verified_per_arm"),
            "paired_to_normal_rho": receipt.get("paired_to_normal_rho"),
        }
        if args.second_host is not None:
            replay_src = args.second_host / f"{case}.json"
            if replay_src.exists():
                if receipt["status"] == "PASS":
                    assert equivalent(json.loads(replay_src.read_text()), receipt), case
                replay_dst = args.out / "second_host_replay" / f"{case}.json"
                shutil.copyfile(replay_src, replay_dst)
                cases[case]["second_host_replay_path"] = str(replay_dst.relative_to(args.out))
                cases[case]["second_host_replay_sha256"] = sha(replay_dst)
    manifest = {
        "schema": "compact-swap-quotient-hosted-archive-v1",
        "source_head": args.head,
        "source_head_kind": "PR head; Actions compiled the recorded PR merge checkout",
        "checkout_merge_commit": args.checkout_merge,
        "checkout_main_parent": args.main_parent,
        "github_run_url": args.run_url,
        "frozen_sha256": sha(args.frozen),
        "cases": cases,
        "extraction": "tar -xzf raw/n41_L1024.tar.gz -C DEST; repeat for n53",
        "independent_replay": "python3 verify_panel.py --n 41 --run-dir DEST/n41_L1024 --out NEW_RECEIPT.json",
    }
    if args.second_host is not None:
        host_src = args.second_host / "HOST.json"
        host_dst = args.out / "second_host_replay" / "HOST.json"
        shutil.copyfile(host_src, host_dst)
        manifest["second_host_replay_host_path"] = str(host_dst.relative_to(args.out))
        manifest["second_host_replay_host_sha256"] = sha(host_dst)
    (args.out / "MANIFEST.json").write_text(json.dumps(
        manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({case: entry["status"] for case, entry in cases.items()},
                     sort_keys=True))


if __name__ == "__main__":
    main()
