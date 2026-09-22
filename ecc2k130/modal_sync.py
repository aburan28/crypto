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

Two things are refused before anything is uploaded, because once uploaded they
are indistinguishable from good points and cost the campaign silently:

  * A run id an AWS slot has used. Seeds are (runId << 48) | (walkIndex << 16)
    and the fleet's slot s runs as run id s + 1, so Modal run r walks AWS slot
    r - 1's trails again and every point it reports is one the store already
    holds under the same seed -- dropped by ON CONFLICT, not a collision, not
    on the page. Runs 1-4 did this against slots 0-3 on 2026-09-19/20. The
    bucket is asked for evidence of the slot (a checkpoint or a dp object) and
    the id is also required to sit in the Modal range. The shared seed
    registry at launch, not this numeric convention, prevents concurrent reuse.
  * A run whose iterations-per-point says it is not walking at the campaign's
    cutoff. Walks stop at their own distinguished point, so a weight-34 walk
    meeting a weight-32 walk is recorded about 12% of the time and a
    weight-35 walk's about 4%; the ratio of the checkpoint header's iterations
    to the corpus's records gives the cutoff away without reading a single
    point. Runs 4242-4245 were at 34 and 35 against the campaign's 32.

Both refusals can be overridden by flag for a run the operator has judged,
and every refusal is logged on every pass so it cannot go quiet.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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

# The campaign these uploads join. Mirrors modal_app.py's constants; the two
# files are deployed separately, so each carries its own copy.
CAMPAIGN_CURVE = 131
MODAL_RUN_ID_MIN = 8000
MODAL_RUN_ID_MAX = 9999
CAMPAIGN_DP_WEIGHT = 32
# Measured on the live fleet at weight 32 (aws/README.md, benchmarks/dp-interval),
# every GPU family agreeing to 0.04 in the exponent. Weights 31 and 33 sit 1.6
# away in either direction, so a tolerance of 1.0 separates neighbours.
CAMPAIGN_ITER_PER_DP_LOG2 = 28.41
DP_RATIO_TOLERANCE_LOG2 = 1.0
# Below this the ratio is dominated by the skew between the header and the
# corpus file (they are fetched seconds apart) rather than by the cutoff.
DP_RATIO_MIN_RECORDS = 50000
# The bucket is asked again for a slot that showed no evidence; one that did
# is remembered for the life of the process.
AWS_EVIDENCE_TTL_S = 600


def log(msg):
    sys.stderr.write(time.strftime("%Y-%m-%dT%H:%M:%SZ ", time.gmtime()) + msg + "\n")
    sys.stderr.flush()


class Policy(object):
    """Which refusals the operator has explicitly lifted for this process."""

    def __init__(self, allow_aws_seed_overlap=False, allow_legacy_run_ids=False,
                 allow_dp_weight_mismatch=False):
        self.allow_aws_seed_overlap = bool(allow_aws_seed_overlap)
        self.allow_legacy_run_ids = bool(allow_legacy_run_ids)
        self.allow_dp_weight_mismatch = bool(allow_dp_weight_mismatch)


DEFAULT_POLICY = Policy()


def aws_slot_prefixes(run_id):
    """Where the AWS fleet leaves traces of the slot that shares this run id."""
    slot = int(run_id) - 1
    if slot < 0:
        return []
    return ["ckpt/slot-%05d" % slot, "ckpt/retired/slot-%05d" % slot, "dp/slot-%05d/" % slot]


def aws_slot_evidence(s3, bucket, run_id):
    """One object key proving an AWS slot has used this run id, or None."""
    for prefix in aws_slot_prefixes(run_id):
        page = s3.list_objects_v2(Bucket=bucket, Prefix=prefix, MaxKeys=1)
        for item in page.get("Contents", []) or []:
            return item["Key"]
    return None


_aws_evidence = {}


def cached_aws_slot_evidence(s3, bucket, run_id, now=None):
    now = time.time() if now is None else now
    hit = _aws_evidence.get(int(run_id))
    if hit is not None:
        key, at = hit
        if key is not None or now - at < AWS_EVIDENCE_TTL_S:
            return key
    key = aws_slot_evidence(s3, bucket, run_id)
    _aws_evidence[int(run_id)] = (key, now)
    return key


def admit_run(s3, bucket, curve, run_id, policy=DEFAULT_POLICY):
    """None when this run may be uploaded into the campaign; else the reason.

    Only campaign-curve runs are judged: other curves are not the campaign's
    corpus. With no S3 client (a dry run) the bucket cannot be asked, and the
    range rule stands on its own.
    """
    if int(curve) != CAMPAIGN_CURVE:
        return None
    run_id = int(run_id)
    if s3 is not None:
        key = cached_aws_slot_evidence(s3, bucket, run_id)
        if key is not None and not policy.allow_aws_seed_overlap:
            return ("run id %d is AWS slot %d's (evidence: s3://%s/%s); the run walks that "
                    "slot's seeds again and every point it reports is a re-report the store "
                    "drops. Stop the run and relaunch it with a run id in %d-%d, or pass "
                    "--allow-aws-seed-overlap to upload anyway."
                    % (run_id, run_id - 1, bucket, key, MODAL_RUN_ID_MIN, MODAL_RUN_ID_MAX))
    if not (MODAL_RUN_ID_MIN <= run_id <= MODAL_RUN_ID_MAX) and not policy.allow_legacy_run_ids:
        return ("run id %d is outside the Modal campaign range %d-%d, where an AWS slot "
                "can claim it; relaunch in range, or pass --allow-legacy-run-ids."
                % (run_id, MODAL_RUN_ID_MIN, MODAL_RUN_ID_MAX))
    return None


def iterations_from_header(header):
    """Total iterations a checkpoint header stands for: per-walk steps x walks."""
    return (int(header["iterBase"]) * int(header["threads"]) * int(header["batch"])
            * max(1, int(header.get("lanes") or 1)))


def estimate_dp_weight(log2_iter_per_dp, m=131):
    """The Hamming-weight cutoff that gives about this many iterations per point.

    Uniform 131-bit strings give 2^29.01 per point at weight 32; the curve's
    x-coordinates measure 2^28.41, and the same 0.6 shift holds at weight 34
    (2^25.84 against the measured 2^25.27), so the binomial tail is shifted by
    it before the nearest weight is read off.
    """
    shift = CAMPAIGN_ITER_PER_DP_LOG2 - theoretical_iter_per_dp_log2(CAMPAIGN_DP_WEIGHT, m)
    best, best_gap = None, None
    for k in range(1, m):
        gap = abs(theoretical_iter_per_dp_log2(k, m) + shift - log2_iter_per_dp)
        if best_gap is None or gap < best_gap:
            best, best_gap = k, gap
    return best


def theoretical_iter_per_dp_log2(weight, m=131):
    total = 0
    for k in range(0, int(weight) + 1):
        total += math.comb(m, k)
    return m - math.log2(total)


def dp_weight_verdict(header, records):
    """None when the run's points-per-iteration fits the campaign cutoff, or
    there is too little to judge; else the reason it does not."""
    if header is None or int(records) < DP_RATIO_MIN_RECORDS:
        return None
    iterations = iterations_from_header(header)
    if iterations <= 0:
        return None
    ratio = math.log2(iterations / float(records))
    if abs(ratio - CAMPAIGN_ITER_PER_DP_LOG2) <= DP_RATIO_TOLERANCE_LOG2:
        return None
    return ("run reports one point per 2^%.2f iterations (%d records over %d iterations); "
            "the campaign's weight-%d cutoff gives 2^%.2f, so this run looks like weight %d. "
            "Walks stop at their own distinguished point, so its meetings with campaign "
            "walks would mostly go unrecorded. Stop it and relaunch at --dp-weight %d, or "
            "pass --allow-dp-weight-mismatch to upload anyway."
            % (ratio, records, iterations, CAMPAIGN_DP_WEIGHT, CAMPAIGN_ITER_PER_DP_LOG2,
               estimate_dp_weight(ratio), CAMPAIGN_DP_WEIGHT))
ORBIT_KEY_RE = re.compile(
    r"^dp/(slot-\d+)/([0-9a-f]{32})-(\d+)-([0-9a-f]{64})\.bin$"
)
CORPUS_RE = re.compile(r"^dp/curve(\d+)-run(\d+)\.bin$")


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


def committed_offset(s3, bucket, run_id, corpus):
    """Recover the contiguous uploaded prefix from immutable S3 objects.

    Local state is a cache, not evidence that points were uploaded. A new
    host or a lost /tmp must not resend a growing corpus from byte zero under
    a new hash. Check every object's digest against this volume snapshot so
    reusing a run id with a different corpus cannot silently skip new points.
    No objects or checkpoints are changed during reconciliation.
    """
    prefix = "dp/slot-%05d/%s-" % (slot_for_run(run_id), stream_id(run_id))
    ranges = []
    request = dict(Bucket=bucket, Prefix=prefix, MaxKeys=1000)
    while True:
        page = s3.list_objects_v2(**request)
        for item in page.get("Contents", []) or []:
            key = item["Key"]
            match = ORBIT_KEY_RE.fullmatch(key)
            if not match:
                if key.endswith(".bin"):
                    raise ValueError("unrecognised corpus object: " + key)
                continue
            start, size = int(match[3]), int(item["Size"])
            if start % RECORD_BYTES or size <= 0 or size % RECORD_BYTES:
                raise ValueError("unaligned corpus object: " + key)
            ranges.append((start, start + size, match[4], key))
        if not page.get("IsTruncated"):
            break
        token = page.get("NextContinuationToken")
        if not token or token == request.get("ContinuationToken"):
            raise ValueError("incomplete S3 listing for " + prefix)
        request["ContinuationToken"] = token

    whole = os.path.getsize(corpus) // RECORD_BYTES * RECORD_BYTES
    offset = 0
    with open(corpus, "rb") as src:
        for start, end, digest, key in sorted(ranges):
            if start > offset:
                raise ValueError("gap in uploaded corpus at byte %d; refusing to skip points" % offset)
            if end > whole:
                raise ValueError("volume corpus is behind S3 at byte %d; wait for a fresh snapshot" % end)
            src.seek(start)
            remaining = end - start
            h = hashlib.sha256()
            while remaining:
                chunk = src.read(min(1024 * 1024, remaining))
                if not chunk:
                    raise ValueError("volume corpus changed during reconciliation")
                h.update(chunk)
                remaining -= len(chunk)
            if h.hexdigest() != digest:
                raise ValueError("volume corpus differs from uploaded object " + key)
            offset = max(offset, end)
    return offset


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


FULL_CHECKPOINT_RETRY_S = 300


def fetch_checkpoint(volume, curve, run_id, state, tmp):
    """The freshest checkpoint artefact on the volume: (local, remote, header).

    Prefer the 40-byte `.hdr` sidecar: ingest's checkpointWork only reads the
    first 40 bytes, and hashing/uploading the 370MB walk state every sync pass
    is how a 60s cadence still left the page idle. If the sidecar is missing,
    fall back to the full file at most once per FULL_CHECKPOINT_RETRY_S.
    (None, None, None) when there is nothing to fetch this pass.
    """
    remote_ck = remote_checkpoint(curve, run_id)
    remote_hdr = remote_checkpoint_header(curve, run_id)
    local_hdr = os.path.join(tmp, os.path.basename(remote_hdr))
    local_ck = os.path.join(tmp, os.path.basename(remote_ck))
    if modal_volume_get(volume, remote_hdr, local_hdr):
        local, remote = local_hdr, remote_hdr
    else:
        last_full = float(state.get("full_checkpoint_check") or 0)
        if state.get("checkpoint_sha256") and (time.time() - last_full) < FULL_CHECKPOINT_RETRY_S:
            log("no status header yet at %s; not re-fetching the 370MB checkpoint"
                % remote_hdr)
            return None, None, None
        if not modal_volume_get(volume, remote_ck, local_ck):
            log("no checkpoint yet at %s or %s on volume %s"
                % (remote_hdr, remote_ck, volume))
            return None, None, None
        local, remote = local_ck, remote_ck
        state["full_checkpoint_check"] = time.time()
    return local, remote, read_checkpoint_header(local)


def upload_checkpoint(s3, bucket, curve, run_id, state, local, remote, dry_run=False):
    """Upload the fetched checkpoint when its content hash changed."""
    if local is None:
        return dict(checkpoint_uploaded=False, checkpoint_key=state.get("checkpoint_key"))
    slot = slot_for_run(run_id)
    digest = sha256_file(local)
    if state.get("checkpoint_sha256") == digest:
        log("checkpoint %s unchanged (%d bytes)" % (remote, os.path.getsize(local)))
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
    return dict(checkpoint_uploaded=True, checkpoint_key=key)


def sync_checkpoint(s3, bucket, volume, curve, run_id, state_dir, dry_run=False):
    """Fetch and upload the checkpoint on its own (no corpus, no admission)."""
    state_file = state_path(curve, run_id, state_dir)
    state = load_state(state_file)
    with tempfile.TemporaryDirectory(prefix="ecc-modal-ckpt-") as tmp:
        local, remote, _ = fetch_checkpoint(volume, curve, run_id, state, tmp)
        result = upload_checkpoint(s3, bucket, curve, run_id, state, local, remote, dry_run=dry_run)
    if not dry_run:
        save_state(state_file, state)
    return result


def refused(run_id, offset, reason):
    return dict(refused=reason, run_id=int(run_id), offset=offset, uploaded_records=0,
                objects=0, pending=False, checkpoint_uploaded=False, checkpoint_key=None)


def sync_once(s3, bucket, volume, curve, run_id, state_dir, dry_run=False,
              policy=DEFAULT_POLICY):
    remote = remote_corpus(curve, run_id)
    state_file = state_path(curve, run_id, state_dir)
    state = load_state(state_file)
    offset = int(state.get("offset", 0))

    # Admission first: nothing of a refused run reaches the bucket, not even
    # its checkpoint header, which would count its work on the page.
    reason = admit_run(s3, bucket, curve, run_id, policy)
    if reason is not None:
        log("REFUSED run %d: %s" % (run_id, reason))
        return refused(run_id, offset, reason)

    with tempfile.TemporaryDirectory(prefix="ecc-modal-sync-") as tmp:
        local = os.path.join(tmp, os.path.basename(remote))
        have_corpus = modal_volume_get(volume, remote, local)
        size = os.path.getsize(local) if have_corpus else 0
        whole = size - size % RECORD_BYTES

        ck_local, ck_remote, header = fetch_checkpoint(volume, curve, run_id, state, tmp)
        if int(curve) == CAMPAIGN_CURVE and header is not None:
            verdict = dp_weight_verdict(header, whole // RECORD_BYTES)
            if verdict is not None and not policy.allow_dp_weight_mismatch:
                log("REFUSED run %d: %s" % (run_id, verdict))
                return refused(run_id, offset, verdict)
            if verdict is not None:
                log("WARNING run %d uploaded anyway (--allow-dp-weight-mismatch): %s"
                    % (run_id, verdict))

        if not have_corpus:
            log("no corpus yet at %s on volume %s (offset still %d)" % (remote, volume, offset))
            ckpt = upload_checkpoint(s3, bucket, curve, run_id, state, ck_local, ck_remote, dry_run)
            if not dry_run:
                save_state(state_file, state)
            return dict(offset=offset, uploaded_records=0, objects=0, pending=False, **ckpt)

        if s3 is not None:
            committed = committed_offset(s3, bucket, run_id, local)
            if committed != offset:
                log("run %d: recovered upload offset %d from S3 (local cache was %d)"
                    % (run_id, committed, offset))
            offset = committed
            state.update(offset=offset, bucket=bucket, remote=remote,
                         slot=slot_for_run(run_id), stream_id=stream_id(run_id))
            # Save point progress independently of checkpoint publication. A
            # failed checkpoint upload must not undo a successful DP upload.
            if not dry_run:
                save_state(state_file, state)

        if whole <= offset:
            log("corpus %s: %d bytes on volume, nothing new after offset %d"
                % (remote, size, offset))
            ckpt = upload_checkpoint(s3, bucket, curve, run_id, state, ck_local, ck_remote, dry_run)
            if not dry_run:
                save_state(state_file, state)
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
        # Points before the checkpoint that follows them, as worker.py does, so
        # the work counter never leads the corpus it accounts for.
        ckpt = upload_checkpoint(s3, bucket, curve, run_id, state, ck_local, ck_remote, dry_run)
        if not dry_run:
            save_state(state_file, state)
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


def sync_pass(s3, bucket, volume, curve, run_ids, state_dir, dry_run=False,
              policy=DEFAULT_POLICY):
    """Sync every run id once. Failures are logged and skipped."""
    results = {}
    for run_id in run_ids:
        try:
            results[run_id] = sync_once(
                s3, bucket, volume, curve, run_id, state_dir, dry_run=dry_run, policy=policy)
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
    ap.add_argument("--allow-aws-seed-overlap", action="store_true",
                    help="upload a run whose run id an AWS slot has used (its points "
                         "replay that slot's; only for finishing a judged run)")
    ap.add_argument("--allow-legacy-run-ids", action="store_true",
                    help="upload a campaign run whose run id is outside %d-%d"
                         % (MODAL_RUN_ID_MIN, MODAL_RUN_ID_MAX))
    ap.add_argument("--allow-dp-weight-mismatch", action="store_true",
                    help="upload a run whose iterations-per-point says it is not at "
                         "the campaign's weight-%d cutoff" % CAMPAIGN_DP_WEIGHT)
    args = ap.parse_args(argv)
    policy = Policy(allow_aws_seed_overlap=args.allow_aws_seed_overlap,
                    allow_legacy_run_ids=args.allow_legacy_run_ids,
                    allow_dp_weight_mismatch=args.allow_dp_weight_mismatch)

    bucket = args.bucket or default_bucket()
    s3 = None
    if not args.dry_run:
        import boto3
        s3 = boto3.client("s3")
    else:
        # A dry run still wants the seed-overlap check when it can have it; the
        # listing writes nothing. Without credentials the range rule stands alone.
        try:
            import boto3
            s3 = boto3.client("s3")
            s3.list_objects_v2(Bucket=bucket, MaxKeys=1)
        except Exception as exc:
            log("dry-run: cannot ask the bucket about AWS slots (%s: %s); only the "
                "run-id range is checked" % (type(exc).__name__, exc))
            s3 = None

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
                s3, bucket, args.volume, args.curve, run_ids, args.state_dir,
                dry_run=args.dry_run, policy=policy)
            if args.watch <= 0:
                print(json.dumps(results, indent=2))
                return 1 if any(r.get("refused") or r.get("error")
                                for r in results.values()) else 0
            for run_id, result in results.items():
                if result.get("uploaded_records") or result.get("checkpoint_uploaded"):
                    print(json.dumps({run_id: result}, indent=2), flush=True)
        if args.watch <= 0:
            return 0
        time.sleep(args.watch)


if __name__ == "__main__":
    raise SystemExit(main())
