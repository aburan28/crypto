#!/usr/bin/env python3
"""Pack hosted raw jobs into deterministic Git artifacts and a hash manifest."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import tarfile


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
                    archive.addfile(info, __import__("io").BytesIO(payload))


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
    assert not args.out.exists(), "refusing to overwrite an evidence archive"
    args.out.mkdir(parents=True)
    (args.out / "raw").mkdir()
    (args.out / "receipts").mkdir()
    (args.out / "hosts").mkdir()
    if args.second_host is not None:
        (args.out / "second_host_replay").mkdir()
    cases = {}
    for n in (37, 41, 53):
        for length in (1, 1024):
            case = f"n{n}_L{length}"
            source = args.runs / case
            assert source.is_dir(), case
            raw = args.out / "raw" / f"{case}.tar.gz"
            pack(source, raw)
            receipt_src = source / "receipt.json"
            receipt = json.loads(receipt_src.read_text())
            assert receipt["status"] == "PASS", case
            assert receipt["n"] == n and receipt["L"] == length, case
            host = json.loads((source / "host.json").read_text())
            assert host["git_head"] == args.checkout_merge, case
            receipt_dst = args.out / "receipts" / f"{case}.json"
            host_dst = args.out / "hosts" / f"{case}.json"
            shutil.copyfile(receipt_src, receipt_dst)
            shutil.copyfile(source / "host.json", host_dst)
            cases[case] = {
                "status": receipt["status"],
                "raw_path": str(raw.relative_to(args.out)),
                "raw_bytes": raw.stat().st_size,
                "raw_sha256": sha(raw),
                "receipt_path": str(receipt_dst.relative_to(args.out)),
                "receipt_sha256": sha(receipt_dst),
                "host_path": str(host_dst.relative_to(args.out)),
                "host_sha256": sha(host_dst),
                "k": receipt.get("k"),
                "targets_verified_per_arm": receipt.get("targets_verified_per_arm"),
                "paired_to_normal_rho": receipt.get("paired_to_normal_rho"),
            }
            if args.second_host is not None:
                replay_src = args.second_host / f"{case}.json"
                assert json.loads(replay_src.read_text()) == receipt, case
                replay_dst = args.out / "second_host_replay" / f"{case}.json"
                shutil.copyfile(replay_src, replay_dst)
                cases[case]["second_host_replay_path"] = str(replay_dst.relative_to(args.out))
                cases[case]["second_host_replay_sha256"] = sha(replay_dst)
    manifest = {
        "schema": "compact-orbit-strong-rho-hosted-archive-v1",
        "source_head": args.head,
        "source_head_kind": "PR head; Actions compiled the recorded PR merge checkout",
        "checkout_merge_commit": args.checkout_merge,
        "checkout_main_parent": args.main_parent,
        "github_run_url": args.run_url,
        "frozen_sha256": sha(args.frozen),
        "cases": cases,
        "extraction": "tar -xzf raw/n37_L1.tar.gz -C DEST; repeat per case",
        "independent_replay": "python3 verify_panel.py --n 37 --L 1 --run-dir DEST/n37_L1 --out NEW_RECEIPT.json",
    }
    if args.second_host is not None:
        second_host_src = args.second_host / "HOST.json"
        second_host_dst = args.out / "second_host_replay" / "HOST.json"
        shutil.copyfile(second_host_src, second_host_dst)
        manifest["second_host_replay_host_path"] = str(second_host_dst.relative_to(args.out))
        manifest["second_host_replay_host_sha256"] = sha(second_host_dst)
    (args.out / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({case: entry["status"] for case, entry in cases.items()},
                     sort_keys=True))


if __name__ == "__main__":
    main()
