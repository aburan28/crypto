#!/usr/bin/env python3
"""Upload new distinguished points from the Modal volume into the campaign bucket.

Modal search appends 32-byte records to /data/dp/curve{C}-run{R}.bin and
checkpoints to /data/ckpt/curve{C}-run{R}.ck on the ecc2k130 volume. The
public dashboard's iteration total and walk rate come from slot checkpoints in
s3://<bucket>/ckpt/, not from the point count. This script copies each new
whole-record stretch into the ecc2k-seed-orbit-v1 dp layout the ingester
already understands, and uploads fresh checkpoints as immutable
ckpt/slot-NNNNN/<sha256>.ck objects so dp_ingest.py's checkpointWork sees
Modal as slot MODAL_SLOT_BASE + run_id.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import struct
import subprocess
import sys
import tempfile
import time

RECORD_BYTES = 32
MODAL_SLOT_BASE = 90000
CKPT_MAGIC = b"ECC2K130"
CKPT_HEADER_BYTES = 40
ORBIT_KEY_RE = re.compile(
    r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$"
)
CORPUS_RE = re.compile(r"^dp/curve(\d+)-run(\d+)\.bin$")


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


def remote_checkpoint(curve, run_id):
    return "ckpt/curve%d-run%d.ck" % (curve, run_id)


def remote_checkpoint_header(curve, run_id):
    """40-byte sidecar. Ingest only reads the header; the 370MB walk state stays on the volume."""
    return "ckpt/curve%d-run%d.hdr" % (curve, run_id)


def pack_checkpoint_header(version, m, threads, batch, lanes, run_id, iter_base):
    return struct.pack("<8s6IQ", CKPT_MAGIC, int(version), int(m), int(threads),
                       int(batch), int(lanes), int(run_id), int(iter_base))


def unpack_checkpoint_header(blob):
    if len(blob) < CKPT_HEADER_BYTES or blob[:8] != CKPT_MAGIC:
        return None
    version, m, threads, batch, lanes, run_id = struct.unpack_from("<6I", blob, 8)
    iter_base, = struct.unpack_from("<Q", blob, 32)
    return dict(version=version, m=m, threads=threads, batch=batch,
                lanes=lanes, runId=run_id, iterBase=iter_base)


def read_checkpoint_header(path):
    try:
        with open(path, "rb") as fh:
            return unpack_checkpoint_header(fh.read(CKPT_HEADER_BYTES))
    except OSError:
        return None


def write_checkpoint_header_file(path, header):
    """Atomic 40-byte file. Same layout as the start of a real checkpoint."""
    if isinstance(header, dict):
        blob = pack_checkpoint_header(
            header["version"], header["m"], header["threads"], header["batch"],
            header["lanes"], header["runId"], header["iterBase"])
    else:
        blob = header
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(blob)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return blob


def iter_base_from_progress(resume_iter_base, pass_iters, walks):
    """Client progress counts (iterBase - start) * walks; invert that."""
    if int(walks) <= 0:
        return int(resume_iter_base)
    return int(resume_iter_base) + int(pass_iters) // int(walks)


def checkpoint_key(slot, ckpt_path):
    return "ckpt/slot-%05d/%s.ck" % (slot, sha256_file(ckpt_path))


def orbit_key(slot, run_id, offset, delta_path):
    sid = stream_id(run_id)
    return "dp/slot-%05d/%s-%016d-%s.bin" % (
        slot, sid, offset, sha256_file(delta_path))


def parse_run_ids(text):
    if not text:
        return []
    return sorted({int(part.strip()) for part in text.split(",") if part.strip()})


def discover_run_ids(volume, curve):
    """Run ids with a corpus on the Modal volume for this curve."""
    proc = subprocess.run(
        ["modal", "volume", "ls", volume, "dp/"],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError("modal volume ls failed: " + err)
    run_ids = []
    for line in (proc.stdout or "").splitlines():
        line = line.strip()
        match = CORPUS_RE.match(line)
        if match and int(match.group(1)) == curve:
            run_ids.append(int(match.group(2)))
    return sorted(set(run_ids))


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


def sync_checkpoint(s3, bucket, volume, curve, run_id, state_dir, dry_run=False):
    """Upload the Modal checkpoint when its content hash changes.

    Prefer the 40-byte `.hdr` sidecar: ingest's checkpointWork only reads the
    first 40 bytes, and hashing/uploading the 370MB walk state every sync pass
    is how a 60s cadence still left the page idle. If the sidecar is missing,
    fall back to the full file at most once per FULL_CHECKPOINT_RETRY_S.
    """
    FULL_CHECKPOINT_RETRY_S = 300
    state_file = state_path(curve, run_id, state_dir)
    state = load_state(state_file)
    remote_ck = remote_checkpoint(curve, run_id)
    remote_hdr = remote_checkpoint_header(curve, run_id)
    slot = slot_for_run(run_id)

    with tempfile.TemporaryDirectory(prefix="ecc-modal-ckpt-") as tmp:
        local_hdr = os.path.join(tmp, os.path.basename(remote_hdr))
        local_ck = os.path.join(tmp, os.path.basename(remote_ck))
        if modal_volume_get(volume, remote_hdr, local_hdr):
            local, remote = local_hdr, remote_hdr
        else:
            last_full = float(state.get("full_checkpoint_check") or 0)
            if state.get("checkpoint_sha256") and (time.time() - last_full) < FULL_CHECKPOINT_RETRY_S:
                log("no status header yet at %s; not re-fetching the 370MB checkpoint"
                    % remote_hdr)
                return dict(checkpoint_uploaded=False, checkpoint_key=state.get("checkpoint_key"))
            if not modal_volume_get(volume, remote_ck, local_ck):
                log("no checkpoint yet at %s or %s on volume %s"
                    % (remote_hdr, remote_ck, volume))
                return dict(checkpoint_uploaded=False, checkpoint_key=None)
            local, remote = local_ck, remote_ck
            state["full_checkpoint_check"] = time.time()

        digest = sha256_file(local)
        if state.get("checkpoint_sha256") == digest:
            log("checkpoint %s unchanged (%d bytes)" % (remote, os.path.getsize(local)))
            if remote == remote_ck and not dry_run:
                save_state(state_file, state)
            return dict(checkpoint_uploaded=False, checkpoint_key=state.get("checkpoint_key"))

        key = checkpoint_key(slot, local)
        if dry_run:
            log("dry-run: would upload checkpoint to s3://%s/%s" % (bucket, key))
        else:
            if s3 is None:
                import boto3
                s3 = boto3.client("s3")
            s3.upload_file(local, bucket, key)
            log("uploaded checkpoint (%d bytes) to s3://%s/%s"
                % (os.path.getsize(local), bucket, key))
        state.update(checkpoint_sha256=digest,
                     checkpoint_key=key,
                     checkpoint_sync=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                     checkpoint_slot=slot,
                     checkpoint_remote=remote)
        if not dry_run:
            save_state(state_file, state)
        return dict(checkpoint_uploaded=True, checkpoint_key=key)


def sync_once(s3, bucket, volume, curve, run_id, state_dir, dry_run=False):
    remote = remote_corpus(curve, run_id)
    state_file = state_path(curve, run_id, state_dir)
    state = load_state(state_file)
    offset = int(state.get("offset", 0))

    with tempfile.TemporaryDirectory(prefix="ecc-modal-sync-") as tmp:
        local = os.path.join(tmp, os.path.basename(remote))
        if not modal_volume_get(volume, remote, local):
            log("no corpus yet at %s on volume %s (offset still %d)" % (remote, volume, offset))
            ckpt = sync_checkpoint(s3, bucket, volume, curve, run_id, state_dir, dry_run=dry_run)
            return dict(offset=offset, uploaded_records=0, objects=0, pending=False, **ckpt)

        size = os.path.getsize(local)
        whole = size - size % RECORD_BYTES
        if whole <= offset:
            log("corpus %s: %d bytes on volume, nothing new after offset %d"
                % (remote, size, offset))
            ckpt = sync_checkpoint(s3, bucket, volume, curve, run_id, state_dir, dry_run=dry_run)
            return dict(offset=offset, uploaded_records=0, objects=0, pending=False, **ckpt)

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
        ckpt = sync_checkpoint(s3, bucket, volume, curve, run_id, state_dir, dry_run=dry_run)
        return dict(offset=offset, uploaded_records=uploaded_records, objects=objects,
                    pending=True, **ckpt)


def run_ids_for_pass(args):
    """Which Modal run ids to sync on this pass."""
    if args.all_runs:
        discovered = discover_run_ids(args.volume, args.curve)
        extra = parse_run_ids(args.run_ids)
        return sorted(set(discovered) | set(extra))
    if args.run_ids:
        return parse_run_ids(args.run_ids)
    return [args.run_id]


def sync_pass(s3, bucket, volume, curve, run_ids, state_dir, dry_run=False):
    """Sync every run id once. Failures are logged and skipped."""
    results = {}
    for run_id in run_ids:
        try:
            results[run_id] = sync_once(
                s3, bucket, volume, curve, run_id, state_dir, dry_run=dry_run)
        except Exception as exc:
            log("sync run %d failed: %s: %s" % (run_id, type(exc).__name__, exc))
            results[run_id] = {"error": str(exc), "run_id": run_id}
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--volume", default=os.environ.get("ECC_MODAL_VOLUME", "ecc2k130"))
    ap.add_argument("--bucket", default="", help="defaults to ECC_BUCKET or ecc2k130-<account>")
    ap.add_argument("--curve", type=int, default=int(os.environ.get("CURVE", "131")))
    ap.add_argument("--run-id", type=int, default=int(os.environ.get("RUNID", "1")))
    ap.add_argument("--run-ids", default=os.environ.get("SYNC_RUN_IDS", ""),
                    help="comma-separated run ids; with --all-runs, merged with discovery")
    ap.add_argument("--all-runs", action="store_true",
                    help="discover every curve*-run*.bin on the Modal volume each pass")
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

    while True:
        try:
            run_ids = run_ids_for_pass(args)
        except Exception as exc:
            log("run discovery failed: %s: %s" % (type(exc).__name__, exc))
            run_ids = parse_run_ids(args.run_ids) or [args.run_id]
        if not run_ids:
            log("no run ids to sync for curve %d on volume %s" % (args.curve, args.volume))
        else:
            log("syncing curve %d runs %s from volume %s to s3://%s/dp/"
                % (args.curve, ",".join(str(r) for r in run_ids), args.volume, bucket))
            results = sync_pass(
                s3, bucket, args.volume, args.curve, run_ids, args.state_dir, dry_run=args.dry_run)
            if args.watch <= 0:
                print(json.dumps(results, indent=2))
                return 0
            for run_id, result in results.items():
                if result.get("uploaded_records") or result.get("checkpoint_uploaded"):
                    print(json.dumps({run_id: result}, indent=2), flush=True)
        if args.watch <= 0:
            return 0
        time.sleep(args.watch)


if __name__ == "__main__":
    raise SystemExit(main())
