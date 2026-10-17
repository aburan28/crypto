#!/usr/bin/env python3
"""Package, deploy and enroll workers for the ECC2K-130 upload gateway.

The command intentionally uses boto3 rather than assuming the AWS CLI is
installed.  Enrollment tokens are written to a mode-0600 file and are never
printed to stdout or stored in CloudFormation.
"""

import argparse
import datetime
import hashlib
import io
import json
import os
import re
import secrets
import stat
import sys
import zipfile

import boto3
from botocore.exceptions import ClientError, WaiterError


STACK = "ecc2k130-upload-gateway"
WORKER_RE = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9_.:-]{0,127}$")


def package(root):
    memory = io.BytesIO()
    with zipfile.ZipFile(memory, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in ("__init__.py", "common.py", "broker.py", "validator.py", "server.py"):
            path = os.path.join(root, name)
            info = zipfile.ZipInfo("upload_gateway/" + name,
                                   date_time=(2026, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (stat.S_IFREG | 0o644) << 16
            with open(path, "rb") as source:
                archive.writestr(info, source.read())
    return memory.getvalue()


def stack(cfn):
    try:
        return cfn.describe_stacks(StackName=STACK)["Stacks"][0]
    except ClientError as error:
        if "does not exist" in str(error):
            return None
        raise


def outputs(cfn):
    found = stack(cfn)
    if not found:
        raise SystemExit("stack %s is not deployed" % STACK)
    return {row["OutputKey"]: row["OutputValue"] for row in found.get("Outputs", [])}


def failureEvents(cfn):
    try:
        events = cfn.describe_stack_events(StackName=STACK)["StackEvents"]
    except ClientError:
        return
    for event in events[:20]:
        if event.get("ResourceStatus", "").endswith("FAILED"):
            print("%s %s: %s" % (event.get("LogicalResourceId"),
                                  event.get("ResourceStatus"),
                                  event.get("ResourceStatusReason", "")), file=sys.stderr)


def deploy(args):
    session = boto3.Session(region_name=args.region)
    account = session.client("sts").get_caller_identity()["Account"]
    canonical = args.bucket or "ecc2k130-" + account
    ingress = args.ingress_bucket or "ecc2k130-ingress-" + account
    root = os.path.dirname(os.path.abspath(__file__))
    body = package(root)
    digest = hashlib.sha256(body).hexdigest()
    codeKey = "infra/upload-gateway/%s.zip" % digest
    session.client("s3").put_object(Bucket=canonical, Key=codeKey, Body=body,
                                    ContentType="application/zip",
                                    Metadata={"sha256": digest})
    templateName = "template.yaml" if args.mode == "lambda" else "template-external.yaml"
    with open(os.path.join(root, templateName)) as source:
        template = source.read()
    cfn = session.client("cloudformation")
    parameters = [
        {"ParameterKey": "CanonicalBucketName", "ParameterValue": canonical},
        {"ParameterKey": "IngressBucketName", "ParameterValue": ingress},
        {"ParameterKey": "CodeBucketName", "ParameterValue": canonical},
        {"ParameterKey": "CodeKey", "ParameterValue": codeKey},
    ]
    if args.mode == "external":
        parameters += [
            {"ParameterKey": "StatusBucketName",
             "ParameterValue": args.status_bucket or "ecc2k130-status-" + account},
            {"ParameterKey": "RdsSecurityGroupId",
             "ParameterValue": args.rds_security_group},
        ]
    common = dict(StackName=STACK, TemplateBody=template, Parameters=parameters,
                  Capabilities=["CAPABILITY_NAMED_IAM"],
                  Tags=[{"Key": "Campaign", "Value": "ecc2k130"},
                        {"Key": "ManagedBy", "Value": "cloudformation"}])
    found = stack(cfn)
    if found and found["StackStatus"] == "ROLLBACK_COMPLETE":
        cfn.delete_stack(StackName=STACK)
        cfn.get_waiter("stack_delete_complete").wait(
            StackName=STACK, WaiterConfig={"Delay": 5, "MaxAttempts": 120})
        found = None
    try:
        if found:
            cfn.update_stack(**common)
            action = "update"
        else:
            cfn.create_stack(**common, OnFailure="ROLLBACK")
            action = "create"
    except ClientError as error:
        if "No updates are to be performed" in str(error):
            action = None
        else:
            raise
    if action:
        try:
            cfn.get_waiter("stack_%s_complete" % action).wait(
                StackName=STACK, WaiterConfig={"Delay": 10, "MaxAttempts": 120})
        except WaiterError:
            failureEvents(cfn)
            raise
    out = outputs(cfn)
    result = {"stack": STACK, "mode": args.mode, "codeSha256": digest,
              "ingressBucket": out["IngressBucket"]}
    if "ApiUrl" in out:
        result["apiUrl"] = out["ApiUrl"]
    if "CodeKey" in out:
        result["codeKey"] = out["CodeKey"]
    print(json.dumps(result, indent=2))


def enroll(args):
    if not WORKER_RE.fullmatch(args.worker_id):
        raise SystemExit("worker id must match " + WORKER_RE.pattern)
    session = boto3.Session(region_name=args.region)
    cfn = session.client("cloudformation")
    table = outputs(cfn)["WorkerTokensTable"]
    token = secrets.token_urlsafe(48)
    digest = hashlib.sha256(token.encode()).hexdigest()
    now = int(datetime.datetime.now(datetime.timezone.utc).timestamp())
    expires = now + args.days * 86400
    path = os.path.abspath(args.token_file)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(path, flags, 0o600)
    registered = False
    ddb = session.client("dynamodb")
    try:
        ddb.put_item(
            TableName=table,
            Item={"tokenHash": {"S": digest}, "workerId": {"S": args.worker_id},
                  "enabled": {"BOOL": True}, "createdAt": {"N": str(now)},
                  "expiresAt": {"N": str(expires)}},
            ConditionExpression="attribute_not_exists(tokenHash)")
        registered = True
        with os.fdopen(descriptor, "w") as output:
            descriptor = None
            output.write(token + "\n")
    except Exception:
        if descriptor is not None:
            os.close(descriptor)
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
        if registered:
            ddb.delete_item(TableName=table, Key={"tokenHash": {"S": digest}})
        raise
    print(json.dumps({"workerId": args.worker_id, "tokenFile": path,
                      "expiresAt": expires}, indent=2))


def revoke(args):
    session = boto3.Session(region_name=args.region)
    table = outputs(session.client("cloudformation"))["WorkerTokensTable"]
    with open(args.token_file) as source:
        token = source.read().strip()
    digest = hashlib.sha256(token.encode()).hexdigest()
    session.client("dynamodb").update_item(
        TableName=table, Key={"tokenHash": {"S": digest}},
        UpdateExpression="SET enabled = :no",
        ExpressionAttributeValues={":no": {"BOOL": False}})
    print(json.dumps({"revoked": True, "tokenFile": os.path.abspath(args.token_file)}))


def gatewayCredentials(args):
    session = boto3.Session(region_name=args.region)
    user = outputs(session.client("cloudformation")).get("ExternalGatewayUser")
    if not user:
        raise SystemExit("the deployed stack does not use the external gateway mode")
    iam = session.client("iam")
    active = [row for row in iam.list_access_keys(UserName=user).get("AccessKeyMetadata", [])
              if row.get("Status") == "Active"]
    if active:
        raise SystemExit("external gateway already has an active access key; revoke it before creating another")
    path = os.path.abspath(args.credentials_file)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    key = None
    try:
        key = iam.create_access_key(UserName=user)["AccessKey"]
        with os.fdopen(descriptor, "w") as output:
            descriptor = None
            json.dump({"AWS_ACCESS_KEY_ID": key["AccessKeyId"],
                       "AWS_SECRET_ACCESS_KEY": key["SecretAccessKey"]}, output)
            output.write("\n")
    except Exception:
        if descriptor is not None:
            os.close(descriptor)
        try:
            os.unlink(path)
        except FileNotFoundError:
            pass
        if key is not None:
            iam.delete_access_key(UserName=user, AccessKeyId=key["AccessKeyId"])
        raise
    print(json.dumps({"user": user, "credentialsFile": path}, indent=2))


def info(args):
    session = boto3.Session(region_name=args.region)
    print(json.dumps(outputs(session.client("cloudformation")), indent=2,
                     sort_keys=True))


def parser():
    result = argparse.ArgumentParser()
    result.add_argument("--region", default=os.environ.get("AWS_DEFAULT_REGION", "us-west-2"))
    sub = result.add_subparsers(dest="command", required=True)
    install = sub.add_parser("deploy")
    install.add_argument("--bucket")
    install.add_argument("--ingress-bucket")
    install.add_argument("--status-bucket")
    install.add_argument("--rds-security-group", default="sg-08c045a6f2fd2b8dc")
    install.add_argument("--mode", choices=("lambda", "external"), default="lambda")
    install.set_defaults(func=deploy)
    add = sub.add_parser("enroll")
    add.add_argument("worker_id")
    add.add_argument("--days", type=int, default=30)
    add.add_argument("--token-file", required=True)
    add.set_defaults(func=enroll)
    remove = sub.add_parser("revoke")
    remove.add_argument("--token-file", required=True)
    remove.set_defaults(func=revoke)
    gateway = sub.add_parser("gateway-credentials")
    gateway.add_argument("--credentials-file", required=True)
    gateway.set_defaults(func=gatewayCredentials)
    show = sub.add_parser("info")
    show.set_defaults(func=info)
    return result


if __name__ == "__main__":
    parsed = parser().parse_args()
    parsed.func(parsed)
