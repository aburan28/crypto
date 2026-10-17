"""HTTP API Lambda for authenticated worker reads, leases and upload tickets."""

import base64
import boto3
import hashlib
import hmac
import json
import os
import secrets
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor

from botocore.exceptions import ClientError
from botocore.config import Config

from .common import (CHECKPOINT_KEY_RE, DP_KEY_RE, LEGACY_CHECKPOINT_KEY_RE,
                     SLOT_KEY_RE, ValidationError, campaignIdentity,
                     canonicalJson, objectSpec, validateTicket)


REGION = os.environ.get("AWS_REGION", "us-west-2")
CANONICAL_BUCKET = os.environ["CANONICAL_BUCKET"]
INGRESS_BUCKET = os.environ["INGRESS_BUCKET"]
AUTH_TABLE = os.environ["AUTH_TABLE"]
TICKET_TABLE = os.environ["TICKET_TABLE"]
VALIDATOR_QUEUE_URL = os.environ["VALIDATOR_QUEUE_URL"]
TICKET_SECONDS = int(os.environ.get("TICKET_SECONDS", "600"))
LEASE_SECONDS = 180
MAX_SLOT = 65534

# Botocore keeps legacy SigV2 as the S3 presign default in some regions.  A
# global SigV2 PUT can answer with a redirect that curl considers successful
# even though no object was stored.  Bind tickets to regional SigV4 URLs.
s3 = boto3.client("s3", region_name=REGION,
                  config=Config(signature_version="s3v4",
                                s3={"addressing_style": "virtual"}))
ddb = boto3.client("dynamodb", region_name=REGION)
sqs = boto3.client("sqs", region_name=REGION)


def response(status, body):
    return {"statusCode": status,
            "headers": {"content-type": "application/json",
                        "cache-control": "no-store"},
            "body": json.dumps(body, sort_keys=True, separators=(",", ":"))}


def jsonBody(event):
    raw = event.get("body") or "{}"
    if event.get("isBase64Encoded"):
        raw = base64.b64decode(raw).decode()
    body = json.loads(raw)
    if not isinstance(body, dict):
        raise ValidationError("request body must be an object")
    return body


def tokenRecord(event):
    headers = {str(k).lower(): v for k, v in (event.get("headers") or {}).items()}
    auth = headers.get("authorization") or ""
    if not auth.startswith("Bearer "):
        raise PermissionError("missing bearer token")
    token = auth[7:]
    if len(token) < 32 or len(token) > 256:
        raise PermissionError("invalid bearer token")
    digest = hashlib.sha256(token.encode()).hexdigest()
    item = ddb.get_item(TableName=AUTH_TABLE,
                        Key={"tokenHash": {"S": digest}},
                        ConsistentRead=True).get("Item")
    now = int(time.time())
    if (not item or item.get("enabled", {}).get("BOOL") is not True
            or int(item.get("expiresAt", {}).get("N", "0")) <= now):
        raise PermissionError("unknown or expired bearer token")
    return {"workerId": item["workerId"]["S"], "tokenHash": digest}


def assertWorker(body, token):
    if body.get("workerId") != token["workerId"]:
        raise PermissionError("token is not enrolled for this worker")


def currentCampaign():
    obj = s3.get_object(Bucket=CANONICAL_BUCKET, Key="campaign.json")
    raw = obj["Body"].read()
    config = json.loads(raw)
    return config, campaignIdentity(config)


def slotKey(slot):
    if isinstance(slot, bool) or not isinstance(slot, int) or not 0 <= slot <= MAX_SLOT:
        raise ValidationError("slot is outside the campaign range")
    return "slots/slot-%05d.json" % slot


def getJson(key):
    try:
        obj = s3.get_object(Bucket=CANONICAL_BUCKET, Key=key)
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") in ("NoSuchKey", "404"):
            return None, None
        raise
    return json.loads(obj["Body"].read()), obj["ETag"]


def putJson(key, value, etag=None, create=False):
    args = {"Bucket": CANONICAL_BUCKET, "Key": key,
            "Body": canonicalJson(value), "ContentType": "application/json"}
    if create:
        args["IfNoneMatch"] = "*"
    elif etag:
        args["IfMatch"] = etag
    try:
        s3.put_object(**args)
        return True
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") in (
                "PreconditionFailed", "ConditionalRequestConflict", "412"):
            return False
        raise


def assertLease(slot, owner, campaignId=None):
    item, _ = getJson(slotKey(slot))
    now = int(time.time())
    if not item or item.get("owner") != owner or int(item.get("leaseUntil") or 0) < now:
        raise PermissionError("worker does not hold the active slot lease")
    if campaignId and item.get("campaignId") not in (None, campaignId):
        raise PermissionError("slot belongs to a different campaign")
    return item


def ticketItem(finalKey):
    return ddb.get_item(TableName=TICKET_TABLE,
                        Key={"finalKey": {"S": finalKey}},
                        ConsistentRead=True).get("Item")


def unpackTicket(item):
    if not item:
        return None
    out = {}
    for key, value in item.items():
        if "S" in value:
            out[key] = value["S"]
        elif "N" in value:
            out[key] = int(value["N"])
        elif "BOOL" in value:
            out[key] = value["BOOL"]
    if "manifest" in out:
        out["manifest"] = json.loads(out["manifest"])
    return out


def issueTicket(body, token):
    assertWorker(body, token)
    spec = validateTicket(body)
    now = int(time.time())
    finalKey = body["key"]
    existing = unpackTicket(ticketItem(finalKey))
    if existing and existing.get("sha256") == body["sha256"]:
        if existing.get("status") == "ACCEPTED":
            return response(200, {"status": "ACCEPTED", "key": finalKey})
        if existing.get("workerId") == token["workerId"]:
            # A restart can move this enrolled worker to another lease while
            # a completed upload is still validating.  Let it recover the
            # terminal/in-flight status without granting a fresh upload URL;
            # requiring the old lease here permanently pins its durable spool.
            if existing.get("status") in ("READY", "VALIDATING", "REJECTED"):
                return ticketResponse(existing)
            if (existing.get("owner") == body["owner"]
                    and existing.get("uploadExpiresAt", 0) > now):
                return ticketResponse(existing)
    config, campaignId = currentCampaign()
    lease = assertLease(spec["slot"], body["owner"], campaignId)
    ticketId = secrets.token_hex(24)
    ingressKey = "incoming/%s/%s" % (ticketId, urllib.parse.quote(finalKey, safe=""))
    expiresAt = now + TICKET_SECONDS
    item = {
        "finalKey": {"S": finalKey}, "ticketId": {"S": ticketId},
        "ingressKey": {"S": ingressKey}, "workerId": {"S": token["workerId"]},
        "owner": {"S": body["owner"]}, "kind": {"S": spec["kind"]},
        "slot": {"N": str(spec["slot"])}, "bytes": {"N": str(body["bytes"])},
        "sha256": {"S": body["sha256"]}, "campaignId": {"S": campaignId},
        "leaseUntil": {"N": str(int(lease["leaseUntil"]))},
        "status": {"S": "UPLOADING"}, "createdAt": {"N": str(now)},
        "uploadExpiresAt": {"N": str(expiresAt)},
        "expiresAt": {"N": str(expiresAt + 86400)},
    }
    for key in ("streamId", "offset"):
        if key in spec:
            kind = "S" if isinstance(spec[key], str) else "N"
            item[key] = {kind: str(spec[key])}
    if body.get("manifest") is not None:
        item["manifest"] = {"S": canonicalJson(body["manifest"]).decode()}
    try:
        ddb.put_item(TableName=TICKET_TABLE, Item=item,
                     ConditionExpression=("attribute_not_exists(finalKey) OR "
                                          "uploadExpiresAt < :now OR #s = :rejected"),
                     ExpressionAttributeNames={"#s": "status"},
                     ExpressionAttributeValues={":now": {"N": str(now)},
                                                ":rejected": {"S": "REJECTED"}})
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") != "ConditionalCheckFailedException":
            raise
        old = unpackTicket(ticketItem(finalKey))
        if old and old.get("sha256") == body["sha256"]:
            return ticketResponse(old)
        raise ValidationError("another ticket already owns this final key")
    return ticketResponse(unpackTicket(item))


def ticketResponse(ticket):
    status = ticket.get("status")
    body = {"status": status, "key": ticket["finalKey"],
            "ticketId": ticket.get("ticketId")}
    if status == "UPLOADING":
        checksum = base64.b64encode(bytes.fromhex(ticket["sha256"])).decode()
        params = {"Bucket": INGRESS_BUCKET, "Key": ticket["ingressKey"],
                  "ContentType": "application/octet-stream",
                  "ChecksumSHA256": checksum,
                  "Metadata": {"ticket-id": ticket["ticketId"]}}
        body["uploadUrl"] = s3.generate_presigned_url(
            "put_object", Params=params, ExpiresIn=TICKET_SECONDS,
            HttpMethod="PUT")
        body["uploadHeaders"] = {
            "content-type": "application/octet-stream",
            "x-amz-checksum-sha256": checksum,
            "x-amz-meta-ticket-id": ticket["ticketId"],
        }
    if ticket.get("reason"):
        body["reason"] = ticket["reason"]
    return response(200, body)


def completeTicket(body, token):
    assertWorker(body, token)
    finalKey = body.get("key")
    ticket = unpackTicket(ticketItem(finalKey))
    if (not ticket or ticket.get("ticketId") != body.get("ticketId")
            or ticket.get("workerId") != token["workerId"]):
        raise PermissionError("ticket does not belong to this worker")
    if ticket["status"] in ("READY", "VALIDATING", "ACCEPTED", "REJECTED"):
        return ticketResponse(ticket)
    if (ticket["status"] != "UPLOADING"
            or ticket.get("uploadExpiresAt", 0) <= int(time.time())):
        raise ValidationError("ticket is not uploadable")
    head = s3.head_object(Bucket=INGRESS_BUCKET, Key=ticket["ingressKey"])
    if head["ContentLength"] != ticket["bytes"]:
        raise ValidationError("uploaded object length does not match its ticket")
    if head.get("Metadata", {}).get("ticket-id") != ticket["ticketId"]:
        raise ValidationError("uploaded object metadata does not match its ticket")
    now = int(time.time())
    ddb.update_item(
        TableName=TICKET_TABLE, Key={"finalKey": {"S": finalKey}},
        UpdateExpression="SET #s = :ready, uploadedAt = :now",
        ConditionExpression="#s = :uploading AND ticketId = :ticket",
        ExpressionAttributeNames={"#s": "status"},
        ExpressionAttributeValues={":ready": {"S": "READY"},
                                   ":uploading": {"S": "UPLOADING"},
                                   ":ticket": {"S": ticket["ticketId"]},
                                   ":now": {"N": str(now)}})
    sqs.send_message(QueueUrl=VALIDATOR_QUEUE_URL,
                     MessageBody=json.dumps({"finalKey": finalKey,
                                             "ticketId": ticket["ticketId"]}))
    ticket["status"] = "READY"
    return ticketResponse(ticket)


def ticketStatus(body, token):
    assertWorker(body, token)
    ticket = unpackTicket(ticketItem(body.get("key")))
    if not ticket or ticket.get("workerId") != token["workerId"]:
        raise PermissionError("ticket does not belong to this worker")
    if body.get("ticketId") and not hmac.compare_digest(
            ticket.get("ticketId", ""), body["ticketId"]):
        raise PermissionError("ticket id does not match")
    return ticketResponse(ticket)


def allowedRead(key):
    if key in ("campaign.json", "solution.json"):
        return True
    if key.startswith("bin/") and ".." not in key:
        return True
    if key in ("aws/worker.py", "aws/protocol.py", "aws/seed_registry.py"):
        return True
    if key.startswith("aws/seed-guard/") and key.endswith("/seed_registry.py") and ".." not in key:
        return True
    payload = key[:-5] if key.endswith(".json") else key
    return (CHECKPOINT_KEY_RE.fullmatch(payload) is not None
            or LEGACY_CHECKPOINT_KEY_RE.fullmatch(payload) is not None
            or DP_KEY_RE.fullmatch(payload) is not None)


def objectAccess(body, token):
    assertWorker(body, token)
    key = body.get("key")
    if not isinstance(key, str) or not allowedRead(key):
        raise PermissionError("object key is outside the worker read allowlist")
    payload = key[:-5] if key.endswith(".json") else key
    if CHECKPOINT_KEY_RE.fullmatch(payload) or LEGACY_CHECKPOINT_KEY_RE.fullmatch(payload):
        spec = objectSpec(payload)
        assertLease(spec["slot"], body.get("owner"))
    if DP_KEY_RE.fullmatch(payload):
        ticket = unpackTicket(ticketItem(payload))
        if not ticket or ticket.get("workerId") != token["workerId"]:
            raise PermissionError("worker may only inspect its own DP uploads")
    try:
        head = s3.head_object(Bucket=CANONICAL_BUCKET, Key=key)
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") in ("NoSuchKey", "404"):
            return response(200, {"exists": False})
        raise
    result = {"exists": True, "bytes": head["ContentLength"],
              "etag": head.get("ETag")}
    if body.get("download"):
        result["downloadUrl"] = s3.generate_presigned_url(
            "get_object", Params={"Bucket": CANONICAL_BUCKET, "Key": key},
            ExpiresIn=300, HttpMethod="GET")
    return response(200, result)


def scanSlots():
    keys = []
    token = None
    while True:
        args = {"Bucket": CANONICAL_BUCKET, "Prefix": "slots/"}
        if token:
            args["ContinuationToken"] = token
        page = s3.list_objects_v2(**args)
        for row in page.get("Contents", []):
            match = SLOT_KEY_RE.fullmatch(row["Key"])
            if not match:
                continue
            keys.append((row["Key"], int(match.group(1))))
        if not page.get("IsTruncated"):
            break
        token = page["NextContinuationToken"]

    def load(row):
        key, slot = row
        item, etag = getJson(key)
        if item is not None:
            item["slot"] = slot
            item["_etag"] = etag
        return item

    # The live campaign has hundreds of small slot documents. Serial S3 GETs
    # exceed a worker's control-plane deadline even though no single request is
    # slow; bounded fan-out keeps claim/scan comfortably below that deadline.
    with ThreadPoolExecutor(max_workers=min(32, max(1, len(keys)))) as pool:
        return [item for item in pool.map(load, keys) if item is not None]


def slotCompatible(item, info, newOnly):
    if item.get("state") not in (None, "active", "idle"):
        return False
    if int(item.get("leaseUntil") or 0) >= int(time.time()):
        return False
    if newOnly and item.get("state") != "idle":
        return False
    family = item.get("gpuFamily")
    return not family or family == info.get("gpuFamily")


def slotClaim(body, token):
    assertWorker(body, token)
    owner = body.get("owner")
    info = body.get("info")
    newOnly = bool(body.get("newOnly"))
    if not isinstance(owner, str) or not owner or not isinstance(info, dict):
        raise ValidationError("slot claim lacks owner or device information")
    campaignId = currentCampaign()[1]
    items = scanSlots()
    for item in sorted(items, key=lambda row: row["slot"]):
        if not slotCompatible(item, info, newOnly):
            continue
        slot, etag = item["slot"], item["_etag"]
        candidate = {key: value for key, value in item.items()
                     if key not in ("slot", "_etag")}
        now = int(time.time())
        candidate.update(owner=owner, leaseUntil=now + LEASE_SECONDS,
                         claimedAt=now, state="active", campaignId=campaignId,
                         workerId=token["workerId"], **info)
        if putJson(slotKey(slot), candidate, etag=etag):
            return response(200, {"slot": slot})
    nextSlot = max((item["slot"] for item in items), default=-1) + 1
    for slot in range(nextSlot, min(nextSlot + 64, MAX_SLOT + 1)):
        now = int(time.time())
        item = dict(owner=owner, leaseUntil=now + LEASE_SECONDS,
                    claimedAt=now, createdAt=now, state="active",
                    campaignId=campaignId, workerId=token["workerId"], **info)
        if putJson(slotKey(slot), item, create=True):
            return response(200, {"slot": slot})
    raise RuntimeError("could not allocate a slot")


def slotHeartbeat(body, token):
    assertWorker(body, token)
    slot, owner = body.get("slot"), body.get("owner")
    fields = body.get("fields") or {}
    if not isinstance(fields, dict):
        raise ValidationError("heartbeat fields must be an object")
    key = slotKey(slot)
    for _ in range(8):
        item, etag = getJson(key)
        now = int(time.time())
        if (not item or item.get("owner") != owner
                or item.get("workerId") not in (None, token["workerId"])
                or int(item.get("leaseUntil") or 0) < now):
            return response(200, {"ok": False})
        item.update(fields, leaseUntil=now + LEASE_SECONDS, updatedAt=now,
                    workerId=token["workerId"])
        if putJson(key, item, etag=etag):
            return response(200, {"ok": True})
    raise RuntimeError("slot heartbeat contention")


def slotRelease(body, token):
    assertWorker(body, token)
    key = slotKey(body.get("slot"))
    for _ in range(8):
        item, etag = getJson(key)
        if (not item or item.get("owner") != body.get("owner")
                or item.get("workerId") not in (None, token["workerId"])):
            return response(200, {"ok": False})
        extra = body.get("extra") or {}
        if not isinstance(extra, dict):
            raise ValidationError("release fields must be an object")
        item.update(extra, leaseUntil=0, state=body.get("state") or "idle")
        if putJson(key, item, etag=etag):
            return response(200, {"ok": True})
    raise RuntimeError("slot release contention")


def seedRead(body, token):
    assertWorker(body, token)
    key = body.get("key")
    if not isinstance(key, str) or not key.startswith("seed-registry/curve-131/") or ".." in key:
        raise PermissionError("seed key is outside the registry")
    value, etag = getJson(key)
    return response(200, {"value": value, "etag": etag})


def seedCas(body, token):
    assertWorker(body, token)
    key = body.get("key")
    value = body.get("value")
    if (not isinstance(key, str) or not key.startswith("seed-registry/curve-131/")
            or ".." in key or not isinstance(value, dict)):
        raise PermissionError("seed update is outside the registry")
    ok = putJson(key, value, etag=body.get("etag"), create=not body.get("etag"))
    return response(200, {"ok": ok})


ROUTES = {
    ("POST", "/v1/uploads/ticket"): issueTicket,
    ("POST", "/v1/uploads/complete"): completeTicket,
    ("POST", "/v1/uploads/status"): ticketStatus,
    ("POST", "/v1/objects/access"): objectAccess,
    ("POST", "/v1/slots/claim"): slotClaim,
    ("POST", "/v1/slots/heartbeat"): slotHeartbeat,
    ("POST", "/v1/slots/release"): slotRelease,
    ("POST", "/v1/slots/scan"): lambda body, token: response(200, {"items": scanSlots()}),
    ("POST", "/v1/seeds/read"): seedRead,
    ("POST", "/v1/seeds/cas"): seedCas,
}


def lambdaHandler(event, context):
    try:
        method = ((event.get("requestContext") or {}).get("http") or {}).get("method")
        method = method or event.get("httpMethod") or ""
        path = event.get("rawPath") or event.get("path") or ""
        handler = ROUTES.get((method.upper(), path))
        if not handler:
            return response(404, {"error": "route not found"})
        token = tokenRecord(event)
        return handler(jsonBody(event), token)
    except PermissionError as error:
        return response(403, {"error": str(error)})
    except (ValidationError, ValueError, json.JSONDecodeError) as error:
        return response(400, {"error": str(error)})
    except ClientError as error:
        code = error.response.get("Error", {}).get("Code", "AWS error")
        return response(503, {"error": code})
    except Exception:
        print("broker request failed", flush=True)
        return response(500, {"error": "internal broker error"})
