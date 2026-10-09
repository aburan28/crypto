"""Archive all frozen paired runs with deterministic, digest-checked receipts."""

import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile


HERE = Path(__file__).resolve().parent
RUN = Path("/private/tmp/s3-c1800-stripe-screen-20261009")
LOCAL = HERE / "local"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    raw_ledger = (RUN / "raw.jsonl").read_bytes()
    rows = [json.loads(line) for line in raw_ledger.splitlines()]
    assert len(rows) == 16
    names = set()
    archive_buffer = io.BytesIO()
    with tarfile.open(fileobj=archive_buffer, mode="w") as archive:
        for row in sorted(rows, key=lambda item: (item["label"], item["arm"], item["repetition"])):
            assert row["verified_native_and_fixture"] is True
            assert row["exit_code"] == 0 and row["error"] is None
            name = f"{row['label']}_{row['arm']}_R{row['repetition']}.json"
            assert name not in names
            names.add(name)
            raw = (RUN / name).read_bytes()
            assert sha(raw) == row["raw_sha256"]
            assert json.loads(raw) == row["report"]
            info = tarfile.TarInfo(name)
            info.size = len(raw)
            info.mode = 0o644
            info.mtime = 0
            archive.addfile(info, io.BytesIO(raw))
    assert len(names) == 16
    LOCAL.mkdir(parents=True, exist_ok=True)
    (LOCAL / "raw_paired16.jsonl.gz").write_bytes(gzip.compress(raw_ledger, mtime=0))
    (LOCAL / "raw_reports.tar.gz").write_bytes(gzip.compress(archive_buffer.getvalue(), mtime=0))
    host = json.loads((RUN / "host.json").read_bytes())
    assert host["host_exclusive_partition_receipt"] is None
    (LOCAL / "host.json").write_text(json.dumps(host, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "rows": len(rows),
        "ledger_decompressed_sha256": sha(raw_ledger),
        "archive_sha256": sha((LOCAL / "raw_reports.tar.gz").read_bytes()),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
