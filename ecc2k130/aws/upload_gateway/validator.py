"""SQS Lambda that admits validated quarantine objects to the campaign corpus."""

import hashlib
import json
import os
import tempfile
import time

import boto3
from botocore.exceptions import ClientError

from .common import (ValidationError, campaignIdentity, canonicalJson,
                     objectSpec, reserveStreamRange, validateDpFile,
                     validateManifest)


REGION = os.environ.get("AWS_REGION", "us-west-2")
CANONICAL_BUCKET = os.environ["CANONICAL_BUCKET"]
INGRESS_BUCKET = os.environ["INGRESS_BUCKET"]
TICKET_TABLE = os.environ["TICKET_TABLE"]
STREAM_TABLE = os.environ["STREAM_TABLE"]

s3 = boto3.client("s3", region_name=REGION)
ddb = boto3.client("dynamodb", region_name=REGION)


def unpack(item):
    out = {}
    for key, value in (item or {}).items():
        if "S" in value:
            out[key] = value["S"]
        elif "N" in value:
            out[key] = int(value["N"])
        elif "BOOL" in value:
            out[key] = value["BOOL"]
    if "manifest" in out:
        out["manifest"] = json.loads(out["manifest"])
    return out


def getTicket(finalKey):
    item = ddb.get_item(TableName=TICKET_TABLE,
                        Key={"finalKey": {"S": finalKey}},
                        ConsistentRead=True).get("Item")
    return unpack(item)


def setStatus(ticket, status, reason=None):
    names = {"#s": "status"}
    values = {":status": {"S": status}, ":now": {"N": str(int(time.time()))},
              ":ticket": {"S": ticket["ticketId"]}}
    expression = "SET #s = :status, validatedAt = :now"
    if reason:
        expression += ", reason = :reason"
        values[":reason"] = {"S": reason[:400]}
    ddb.update_item(TableName=TICKET_TABLE,
                    Key={"finalKey": {"S": ticket["finalKey"]}},
                    UpdateExpression=expression,
                    ConditionExpression="ticketId = :ticket",
                    ExpressionAttributeNames=names,
                    ExpressionAttributeValues=values)


def claimTicket(ticket):
    try:
        ddb.update_item(
            TableName=TICKET_TABLE,
            Key={"finalKey": {"S": ticket["finalKey"]}},
            UpdateExpression="SET #s = :validating, validationStartedAt = :now",
            ConditionExpression="#s = :ready AND ticketId = :ticket",
            ExpressionAttributeNames={"#s": "status"},
            ExpressionAttributeValues={":validating": {"S": "VALIDATING"},
                                       ":ready": {"S": "READY"},
                                       ":ticket": {"S": ticket["ticketId"]},
                                       ":now": {"N": str(int(time.time()))}})
        return True
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") == "ConditionalCheckFailedException":
            return False
        raise


def currentCampaign():
    raw = s3.get_object(Bucket=CANONICAL_BUCKET, Key="campaign.json")["Body"].read()
    config = json.loads(raw)
    return config, campaignIdentity(config)


def validateLease(ticket):
    key = "slots/slot-%05d.json" % ticket["slot"]
    item = json.loads(s3.get_object(Bucket=CANONICAL_BUCKET, Key=key)["Body"].read())
    if item.get("owner") != ticket["owner"]:
        raise ValidationError("slot owner changed before validation")
    if item.get("workerId") not in (None, ticket["workerId"]):
        raise ValidationError("slot belongs to another enrolled worker")
    if int(ticket.get("uploadedAt", 0)) > int(ticket.get("leaseUntil", 0)):
        raise ValidationError("upload completed after the authorized lease")
    if item.get("campaignId") not in (None, ticket["campaignId"]):
        raise ValidationError("slot belongs to another campaign")


def downloadAndHash(ticket, path):
    digest = hashlib.sha256()
    size = 0
    obj = s3.get_object(Bucket=INGRESS_BUCKET, Key=ticket["ingressKey"])
    with open(path, "wb") as output:
        for chunk in iter(lambda: obj["Body"].read(1024 * 1024), b""):
            size += len(chunk)
            if size > ticket["bytes"]:
                raise ValidationError("quarantine object exceeds its ticketed size")
            digest.update(chunk)
            output.write(chunk)
    if size != ticket["bytes"]:
        raise ValidationError("quarantine object size does not match its ticket")
    if digest.hexdigest() != ticket["sha256"]:
        raise ValidationError("quarantine object checksum does not match its ticket")


def reserveOffset(ticket, recordBytes, headerBytes=0):
    if ticket["kind"] != "dp":
        return
    start = int(ticket["offset"])
    streamKey = "%05d:%s" % (ticket["slot"], ticket["streamId"])
    key = {"slotStream": {"S": streamKey}}
    old = unpack(ddb.get_item(TableName=STREAM_TABLE, Key=key,
                              ConsistentRead=True).get("Item"))
    end = start + int(recordBytes)
    if old and old.get("lastKey") == ticket["finalKey"] and old.get("nextOffset") == end:
        return
    end = reserveStreamRange(start, recordBytes, headerBytes,
                             old.get("nextOffset") if old else None)
    values = {":end": {"N": str(end)}, ":key": {"S": ticket["finalKey"]},
              ":now": {"N": str(int(time.time()))}}
    if old:
        values[":start"] = {"N": str(start)}
        condition = "nextOffset = :start"
    else:
        # Framed v2/table-v3 source files begin their record stream after the
        # 16-byte header.  Each uploaded delta repeats that header, but the key
        # offset and nextOffset both describe positions in the source file, so
        # only the first record offset (not zero) is a valid stream start.
        condition = "attribute_not_exists(slotStream)"
    try:
        ddb.update_item(TableName=STREAM_TABLE, Key=key,
                        UpdateExpression="SET nextOffset = :end, lastKey = :key, updatedAt = :now",
                        ConditionExpression=condition,
                        ExpressionAttributeValues=values)
    except ClientError as error:
        if error.response.get("Error", {}).get("Code") == "ConditionalCheckFailedException":
            again = unpack(ddb.get_item(TableName=STREAM_TABLE, Key=key,
                                        ConsistentRead=True).get("Item"))
            if again.get("lastKey") == ticket["finalKey"] and again.get("nextOffset") == end:
                return
            raise ValidationError("DP stream offset was claimed concurrently")
        raise


def publish(ticket, path, validation):
    metadata = {"validated": "true", "sha256": ticket["sha256"],
                "campaign-id": ticket["campaignId"],
                "worker-id": ticket["workerId"], "ticket-id": ticket["ticketId"]}
    with open(path, "rb") as source:
        s3.put_object(Bucket=CANONICAL_BUCKET, Key=ticket["finalKey"],
                      Body=source, ContentType="application/octet-stream",
                      Metadata=metadata)
    manifest = ticket.get("manifest")
    if manifest is not None:
        admitted = dict(manifest)
        admitted.update(validated=True, validatedAt=int(time.time()),
                        ticketId=ticket["ticketId"], records=validation.get("records"))
        s3.put_object(Bucket=CANONICAL_BUCKET, Key=ticket["finalKey"] + ".json",
                      Body=canonicalJson(admitted), ContentType="application/json",
                      Metadata=metadata)


def validateOne(finalKey, ticketId):
    ticket = getTicket(finalKey)
    if not ticket or ticket.get("ticketId") != ticketId:
        return
    if ticket.get("status") == "ACCEPTED":
        return
    if ticket.get("status") == "READY":
        if not claimTicket(ticket):
            return
        ticket = getTicket(finalKey)
    elif ticket.get("status") != "VALIDATING":
        return
    try:
        spec = objectSpec(ticket["finalKey"])
        if spec["kind"] != ticket["kind"] or spec["slot"] != ticket["slot"]:
            raise ValidationError("ticket metadata does not match its final key")
        config, campaignId = currentCampaign()
        if campaignId != ticket["campaignId"]:
            raise ValidationError("campaign changed before validation")
        validateLease(ticket)
        with tempfile.TemporaryDirectory(prefix="ecc-ingress-") as tmp:
            path = os.path.join(tmp, "body")
            downloadAndHash(ticket, path)
            validation = {}
            if ticket["kind"] == "dp":
                validation = validateDpFile(path, ticket["slot"], config)
            elif ticket["kind"] == "solution":
                try:
                    with open(path) as source:
                        value = json.load(source)
                except (OSError, ValueError) as error:
                    raise ValidationError("solution is not valid JSON") from error
                if not isinstance(value, dict) or value.get("slot") not in (None, ticket["slot"]):
                    raise ValidationError("solution does not match its slot")
            validateManifest(ticket.get("manifest"), ticket, campaignId,
                             validation.get("records"))
            reserveOffset(ticket, validation.get("recordBytes", 0),
                          validation.get("headerBytes", 0))
            publish(ticket, path, validation)
        setStatus(ticket, "ACCEPTED")
        s3.delete_object(Bucket=INGRESS_BUCKET, Key=ticket["ingressKey"])
    except ValidationError as error:
        setStatus(ticket, "REJECTED", str(error))


def messages(event):
    for record in event.get("Records", []):
        body = json.loads(record.get("body") or "{}")
        if "finalKey" in body and "ticketId" in body:
            yield body["finalKey"], body["ticketId"]
            continue
        # Native S3 notifications are a recovery signal only.  The explicit
        # completion message names the ticket and is authoritative.
        for item in body.get("Records", []):
            key = item.get("s3", {}).get("object", {}).get("key", "")
            key = key.replace("+", " ")
            if not key.startswith("incoming/"):
                continue
            ticketId = key.split("/", 2)[1]
            # An S3 event can race the completion call.  It has no final key;
            # leave it for the explicit READY message rather than scanning.
            print("quarantine object arrived for ticket %s" % ticketId, flush=True)


def lambdaHandler(event, context):
    for finalKey, ticketId in messages(event):
        validateOne(finalKey, ticketId)
    return {"ok": True}
