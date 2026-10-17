#!/usr/bin/env python3
"""Preserve both finalist cloud artifacts with identical base dumps stored once."""

import argparse
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN_ID = 39005111289
HEAD_SHA = "416489e09c3d8885b6c12f262935e6ba29373460"
SOURCE_SHA = "8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12"
PROVIDERS = (
    (0, 12095116702, "50db2c8e06c884ad9f6a2731aa350dc8bde708d8ea0c52fddf740404874908a4"),
    (1, 12094374158, "737df73989d1a2e663ef3ea9253bb088b6aa92d2734070cae3a3d4b641b75c20"),
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def write_once(path, data):
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError(f"sealed artifact differs: {path}")
    else:
        path.write_bytes(data)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact_root", type=Path)
    args = parser.parse_args()
    root = args.artifact_root.resolve()
    paths = {}
    providers = []
    source_receipts = []
    for half, artifact_id, expected_zip_sha in PROVIDERS:
        prefix = f"n83-e2e-final-half-{half}-{RUN_ID}-1"
        directory = root / prefix
        archive = root / f"provider_half{half}.zip"
        if not directory.is_dir() or sha256(archive) != expected_zip_sha:
            raise ValueError(f"half {half}: provider ZIP or extraction differs")
        analysis = json.loads((directory / "results/analysis.json").read_text())
        replay = json.loads((directory / "local_replay.stdout.txt").read_text())
        source = json.loads((directory / "results/source.json").read_text())
        host = json.loads((directory / "results/host.json").read_text())
        if (analysis != replay or analysis["status"] != "HALF_ADMITTED"
                or source["source_sha256"] != SOURCE_SHA
                or host["git_sha"] != HEAD_SHA
                or (directory / "local_replay.exit_status.txt").read_text().strip() != "0"
                or (directory / "local_replay.stderr.txt").read_bytes()):
            raise ValueError(f"half {half}: independent replay is not admitted")
        source_receipts.append({"source_sha256": source["source_sha256"],
                                "binary_sha256": source["binary_sha256"],
                                "rho_binary_sha256": source["rho_binary_sha256"],
                                "instance_id": host["instance_id"]})
        providers.append({"half": half, "artifact_id": artifact_id,
                          "zip_sha256": expected_zip_sha,
                          "zip_bytes": archive.stat().st_size,
                          "url": ("https://api.github.com/repos/aburan28/crypto/"
                                  f"actions/artifacts/{artifact_id}/zip")})
        for path in sorted(directory.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"unexpected provider symlink: {path}")
            if path.is_file():
                name = f"half{half}/{path.relative_to(directory).as_posix()}"
                paths[name] = path

    entries = {}
    first_base_by_sha = {}
    if (len({item["binary_sha256"] for item in source_receipts}) != 1
            or len({item["rho_binary_sha256"] for item in source_receipts}) != 1
            or len({item["instance_id"] for item in source_receipts}) != 2):
        raise ValueError("host or executed binary strata differ")
    for name, path in sorted(paths.items()):
        digest = sha256(path)
        blob = name
        if name.endswith("/base_dump.jsonl"):
            blob = first_base_by_sha.setdefault(digest, name)
        entries[name] = {"sha256": digest, "bytes": path.stat().st_size,
                         "blob": blob}
    manifest = {"schema": "n83.finalist-dedup-archive/v1",
                "run_id": RUN_ID, "run_attempt": 1, "head_sha": HEAD_SHA,
                "source_sha256": SOURCE_SHA,
                "providers": providers, "source_receipts": source_receipts,
                "entries": entries}
    manifest_bytes = canonical_json(manifest)
    write_once(HERE / "raw_dedup_manifest.json", manifest_bytes)

    output = HERE / "raw_dedup.tar.gz"
    if output.exists():
        raise ValueError(f"archive already exists: {output}")
    temporary = HERE / "raw_dedup.tar.gz.part"
    if temporary.exists():
        raise ValueError(f"incomplete earlier archive exists: {temporary}")
    blobs = {item["blob"] for item in entries.values()}
    with temporary.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0,
                           compresslevel=6) as compressed:
            with tarfile.open(fileobj=compressed, mode="w") as tar:
                for name in sorted(blobs):
                    path = paths[name]
                    info = tarfile.TarInfo(name)
                    info.size = path.stat().st_size
                    info.mode = 0o644
                    info.mtime = 0
                    with path.open("rb") as stream:
                        tar.addfile(info, stream)
                info = tarfile.TarInfo("MANIFEST.json")
                info.size = len(manifest_bytes)
                info.mode = 0o644
                info.mtime = 0
                tar.addfile(info, io.BytesIO(manifest_bytes))
    temporary.replace(output)
    receipt = {"schema": "n83.finalist-dedup-archive-receipt/v1",
               "archive_sha256": sha256(output),
               "archive_bytes": output.stat().st_size,
               "manifest_sha256": sha256(HERE / "raw_dedup_manifest.json"),
               "file_count": len(entries), "stored_blob_count": len(blobs),
               "unique_base_dump_count": len(first_base_by_sha)}
    write_once(HERE / "raw_dedup_receipt.json", canonical_json(receipt))
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
