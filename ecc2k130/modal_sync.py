#!/usr/bin/env python3
"""Upload new distinguished points from the Modal volume into the campaign bucket.

Modal search appends 32-byte records to /data/dp/curve{C}-run{R}.bin on the ecc2k130
volume. The public dashboard reads Postgres, which dp_ingest.py fills from
s3://<bucket>/dp/. This script copies each new whole-record stretch into the
ecc2k-seed-orbit-v1 object layout the ingester already understands, so Modal
contributions ride the same store path as fleet workers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import time

RECORD_BYTES = 32
MODAL_SLOT_BASE = 90000
ORBIT_KEY_RE = re.compile(
    r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$"
)


def log(msg):
    sys.stderr.write(time.strftime("%Y-%m-%dT%H:%M:%SZ ", time.gmtime()) + msg + "\n")
    sys.stderr.flush()


def stream_id(run_id):
    return hashlib.sha256(("modal-run-%d" % run_id).encode()).hexdigest()[:32]


def slot_for_run(run_id):
    slot = MODAL_SLOT_BASE + int(run_id)
    if slot > 99999:
        raise ValueError("run_id %d exceeds modal slot range" % run_id)
    return slot


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def default_bucket():
    bucket = os.environ.get("ECC_BUCKET")
    if bucket:
        return bucket
    import boto3

    account = boto3.client("sts").get_caller_identity()["Account"]
    return "ecc2k130-%s" % account


def state_path(curve, run_id, state_dir):
    os.makedirs(state_dir, exist_ok=True)
    return os.path.join(state_dir, "curve%d-run%d.json" % (curve, run_id))


def load_state(path):
    if not os.path.isfile(path):
        return {"offset": 0, "uploaded_records": 0}
    with open(path) as fh:
        return json.load(fh)


def save_state(path, state):
    tmp = path + ".part"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def remote_corpus(curve, run_id):
    return "dp/curve%d-run%d.bin" % (curve, run_id)


def orbit_key(slot, run_id, offset, delta_path):
    sid = stream_id(run_id)
    return "dp/slot-%05d/%s-%016d-%s.bin" % (
        slot, sid, offset, sha256_file(delta_path))


def modal_volume_get(volume, remote, local):
    r = subprocess.run(
        ["modal", "volume", "get", volume, remote, local],
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        err = (r.stderr or r.stdout or "").strip()
        if "No such file" in err or "not found" in err.lower():
            return False
        raise RuntimeError("modal volume get failed: " + err)
    return True


def sync_once(s3, bucket, volume, curve, run_id, state_dir, dry_run=False):
    remote = remote_corpus(curve, run_id)
    state_file = state_path(curve, run_id, state_dir)
    state = load_state(state_file)
    offset = int(state.get("offset", 0))

    with tempfile.TemporaryDirectory(prefix="ecc-modal-sync-") as tmp:
        local = os.path.join(tmp, os.path.basename(remote))
        if not modal_volume_get(volume, remote, local):
            log("no corpus yet at %s on volume %s (offset still %d)" % (remote, volume, offset))
            return dict(offset=offset, uploaded_records=0, objects=0, pending=False)

        size = os.path.getsize(local)
        whole = size - size % RECORD_BYTES
        if whole <= offset:
            log("corpus %s: %d bytes on volume, nothing new after offset %d"
                % (remote, size, offset))
            return dict(offset=offset, uploaded_records=0, objects=0, pending=False)

        slot = slot_for_run(run_id)
        uploaded_records = 0
        objects = 0
        while offset < whole:
            delta_path = os.path.join(tmp, "delta-%d.bin" % offset)
            with open(local, "rb") as src, open(delta_path, "wb") as out:
                src.seek(offset)
                out.write(src.read(whole - offset))
            key = orbit_key(slot, run_id, offset, delta_path)
            records = os.path.getsize(delta_path) // RECORD_BYTES
            if not ORBIT_KEY_RE.match(key):
                raise RuntimeError("generated key does not match ingest pattern: " + key)
            if dry_run:
                log("dry-run: would upload %d records to s3://%s/%s" % (records, bucket, key))
            else:
                if s3 is None:
                    import boto3
                    s3 = boto3.client("s3")
                s3.upload_file(delta_path, bucket, key)
                log("uploaded %d records to s3://%s/%s" % (records, bucket, key))
            offset += records * RECORD_BYTES
            uploaded_records += records
            objects += 1

        state.update(offset=offset,
                     uploaded_records=int(state.get("uploaded_records", 0)) + uploaded_records,
                     last_sync=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                     slot=slot,
                     stream_id=stream_id(run_id),
                     remote=remote,
                     bucket=bucket)
        if not dry_run:
            save_state(state_file, state)
        return dict(offset=offset, uploaded_records=uploaded_records, objects=objects, pending=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--volume", default=os.environ.get("ECC_MODAL_VOLUME", "ecc2k130"))
    ap.add_argument("--bucket", default="", help="defaults to ECC_BUCKET or ecc2k130-<account>")
    ap.add_argument("--curve", type=int, default=int(os.environ.get("CURVE", "131")))
    ap.add_argument("--run-id", type=int, default=int(os.environ.get("RUNID", "1")))
    ap.add_argument("--state-dir", default=os.environ.get("ECC_MODAL_SYNC_STATE",
                                                          os.path.join(tempfile.gettempdir(),
                                                                       "ecc2k130-modal-sync")))
    ap.add_argument("--watch", type=float, default=0,
                    help="seconds between passes; 0 means one pass")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    bucket = args.bucket or default_bucket()
    s3 = None
    if not args.dry_run:
        import boto3
        s3 = boto3.client("s3")

    log("syncing curve %d run %d from volume %s to s3://%s/dp/"
        % (args.curve, args.run_id, args.volume, bucket))

    while True:
        result = sync_once(s3, bucket, args.volume, args.curve, args.run_id,
                             args.state_dir, dry_run=args.dry_run)
        if args.watch <= 0:
            print(json.dumps(result, indent=2))
            return 0
        if result["uploaded_records"]:
            print(json.dumps(result, indent=2), flush=True)
        time.sleep(args.watch)


if __name__ == "__main__":
    raise SystemExit(main())
