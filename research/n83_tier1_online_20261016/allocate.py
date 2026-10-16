#!/usr/bin/env python3
"""Bounded, receipt-producing AWS Spot allocation for the N83 Tier 1 panel."""

import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import secrets
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNNER = HERE / "run.py"
REGION = "us-east-1"
AZ = "us-east-1b"
INSTANCE_TYPE = "g5g.metal"
SUBNET = "subnet-075981b27daff982c"
VPC = "vpc-0da45a1bf46c94fa9"
MAX_PRICE = "0.50"
MAX_MINUTES = 165
VOLUME_GIB = 50
BRANCH = "codex/n83-tier1-online-20261016"


def command(*args, allow_error=False):
    result = subprocess.run(args, capture_output=True, text=True, check=False,
                            timeout=60)
    if result.returncode and not allow_error:
        raise RuntimeError(f"{args[0]} {args[1]} failed: {result.stderr.strip()}")
    return result


def aws(*args):
    output = command("aws", *args, "--output", "json").stdout.strip()
    return json.loads(output) if output else {}


def aws_try(*args):
    result = command("aws", *args, "--output", "json", allow_error=True)
    if result.returncode:
        return None, result.stderr
    return json.loads(result.stdout) if result.stdout.strip() else {}, None


def current_infra():
    account = aws("sts", "get-caller-identity")["Account"]
    subnet = aws("ec2", "describe-subnets", "--region", REGION,
                 "--subnet-ids", SUBNET)["Subnets"][0]
    if (subnet["VpcId"] != VPC or subnet["AvailabilityZone"] != AZ
            or subnet["State"] != "available"
            or not subnet["MapPublicIpOnLaunch"]):
        raise ValueError("approved subnet or zone differs")
    now = datetime.now(timezone.utc)
    prices = aws("ec2", "describe-spot-price-history", "--region", REGION,
                 "--availability-zone", AZ, "--instance-types", INSTANCE_TYPE,
                 "--product-descriptions", "Linux/UNIX", "--start-time",
                 (now - timedelta(days=1)).isoformat(), "--max-results", "20")
    history = prices["SpotPriceHistory"]
    if not history:
        raise ValueError("current g5g.metal Spot price unavailable")
    history.sort(key=lambda item: item["Timestamp"], reverse=True)
    price = float(history[0]["SpotPrice"])
    if price > float(MAX_PRICE):
        raise ValueError("current Spot price exceeds frozen hourly cap")
    ami = aws("ssm", "get-parameter", "--region", REGION, "--name",
              "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64")[
                  "Parameter"]["Value"]
    image = aws("ec2", "describe-images", "--region", REGION,
                "--image-ids", ami)["Images"][0]
    if image["Architecture"] != "arm64" or image["State"] != "available":
        raise ValueError("Amazon Linux 2023 ARM64 AMI differs")
    return {"account": account, "subnet": SUBNET, "vpc": VPC, "zone": AZ,
            "instance_type": INSTANCE_TYPE, "ami": ami,
            "ami_root_device": image["RootDeviceName"],
            "spot_price": history[0], "checked_utc": now.isoformat(),
            "max_price_per_hour_usd": MAX_PRICE,
            "max_minutes": MAX_MINUTES, "volume_gib": VOLUME_GIB,
            "max_compute_charge_usd": float(MAX_PRICE) * MAX_MINUTES / 60}


def source_preflight(mode):
    spec = importlib.util.spec_from_file_location("n83_tier1_runner", RUNNER)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.fixture_preflight()
    if mode == "selected_panel":
        runner.selection_preflight()
    elif mode == "same_instance_selector":
        if (HERE / "selection.json").exists():
            raise ValueError("host selector allocation requires no prior selection")
        selector_spec = importlib.util.spec_from_file_location(
            "n83_tier1_selector_runner", HERE / "selector_run.py")
        selector = importlib.util.module_from_spec(selector_spec)
        selector_spec.loader.exec_module(selector)
        selector.frozen_inputs()
    else:
        raise ValueError("unsupported allocation mode")
    head = command("git", "-C", str(ROOT), "rev-parse", "HEAD").stdout.strip()
    dirty = command("git", "-C", str(ROOT), "status", "--porcelain",
                    "--untracked-files=all").stdout
    remote = command("git", "-C", str(ROOT), "ls-remote", "origin",
                     f"refs/heads/{BRANCH}").stdout
    if dirty or not remote or remote.split()[0] != head:
        raise ValueError("exact selected source must be clean and pushed before allocation")
    return head


def ensure_role(account):
    name = "n83-tier1-ssm-20261016"
    role, error = aws_try("iam", "get-role", "--role-name", name)
    if role is None:
        if "NoSuchEntity" not in error:
            raise RuntimeError(error)
        trust = {"Version": "2012-10-17", "Statement": [{
            "Effect": "Allow", "Principal": {"Service": "ec2.amazonaws.com"},
            "Action": "sts:AssumeRole"}]}
        aws("iam", "create-role", "--role-name", name,
            "--assume-role-policy-document", json.dumps(trust),
            "--tags", "Key=Project,Value=N83Tier1")
    elif not any(item["Key"] == "Project" and item["Value"] == "N83Tier1"
                 for item in role["Role"].get("Tags", [])):
        raise ValueError("existing IAM role lacks N83 ownership tag")
    aws("iam", "attach-role-policy", "--role-name", name,
        "--policy-arn", "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore")
    bucket = f"aburan28-n83-tier1-20261016-{account}"
    policy = {"Version": "2012-10-17", "Statement": [{
        "Effect": "Allow", "Action": ["s3:PutObject"],
        "Resource": f"arn:aws:s3:::{bucket}/n83-tier1-online-20261016/*"}]}
    aws("iam", "put-role-policy", "--role-name", name,
        "--policy-name", "n83-tier1-artifact-put",
        "--policy-document", json.dumps(policy))
    profile, error = aws_try("iam", "get-instance-profile",
                             "--instance-profile-name", name)
    if profile is None:
        if "NoSuchEntity" not in error:
            raise RuntimeError(error)
        aws("iam", "create-instance-profile", "--instance-profile-name", name,
            "--tags", "Key=Project,Value=N83Tier1")
        aws("iam", "add-role-to-instance-profile", "--instance-profile-name",
            name, "--role-name", name)
    else:
        attached = [item["RoleName"] for item in
                    profile["InstanceProfile"]["Roles"]]
        if not attached:
            aws("iam", "add-role-to-instance-profile", "--instance-profile-name",
                name, "--role-name", name)
        elif attached != [name]:
            raise ValueError("existing instance profile has another role")
    return name, bucket


def ensure_bucket(bucket):
    listed = aws("s3api", "list-buckets")["Buckets"]
    if bucket not in {item["Name"] for item in listed}:
        aws("s3api", "create-bucket", "--region", REGION, "--bucket", bucket)
    aws("s3api", "put-public-access-block", "--bucket", bucket,
        "--public-access-block-configuration", json.dumps({
            "BlockPublicAcls": True, "IgnorePublicAcls": True,
            "BlockPublicPolicy": True, "RestrictPublicBuckets": True}))
    aws("s3api", "put-bucket-encryption", "--bucket", bucket,
        "--server-side-encryption-configuration", json.dumps({
            "Rules": [{"ApplyServerSideEncryptionByDefault": {
                "SSEAlgorithm": "AES256"}}]}))
    aws("s3api", "put-bucket-versioning", "--bucket", bucket,
        "--versioning-configuration", "Status=Enabled")
    if aws("s3api", "get-bucket-versioning", "--bucket", bucket).get(
            "Status") != "Enabled":
        raise ValueError("selection bucket versioning is not enabled")


def ensure_security_group():
    name = "n83-tier1-no-inbound-20261016"
    groups = aws("ec2", "describe-security-groups", "--region", REGION,
                 "--filters", f"Name=vpc-id,Values={VPC}",
                 f"Name=group-name,Values={name}")["SecurityGroups"]
    if groups:
        group = groups[0]
        if group["IpPermissions"]:
            raise ValueError("dedicated N83 group has inbound permissions")
        if not any(permission.get("IpProtocol") == "-1"
                   and any(address.get("CidrIp") == "0.0.0.0/0"
                           for address in permission.get("IpRanges", []))
                   for permission in group["IpPermissionsEgress"]):
            raise ValueError("dedicated N83 group lacks expected outbound route")
        return group["GroupId"]
    created = aws("ec2", "create-security-group", "--region", REGION,
                  "--group-name", name, "--description",
                  "N83 Tier 1 outbound-only research allocation", "--vpc-id", VPC)
    group_id = created["GroupId"]
    aws("ec2", "create-tags", "--region", REGION, "--resources", group_id,
        "--tags", "Key=Project,Value=N83Tier1")
    return group_id


def launch(output, infrastructure, commit, mode):
    output.mkdir(parents=True, exist_ok=False)
    (output / "infrastructure-preflight.json").write_text(json.dumps(
        infrastructure, indent=2, sort_keys=True) + "\n")
    account = infrastructure["account"]
    role, bucket = ensure_role(account)
    ensure_bucket(bucket)
    group = ensure_security_group()
    token = "n83tier1-20261016-" + secrets.token_hex(4)
    deadline = datetime.now(timezone.utc) + timedelta(minutes=MAX_MINUTES)
    template = (HERE / "aws_bootstrap.sh").read_text()
    rendered = (template.replace("__BUCKET__", bucket)
                .replace("__COMMIT__", commit)
                .replace("__RUN_TOKEN__", token)
                .replace("__MODE__", mode))
    user_data = output / "user-data.sh"
    user_data.write_text(rendered)
    command("bash", "-n", str(user_data))
    tags = [{"Key": "Project", "Value": "N83Tier1"},
            {"Key": "ClientToken", "Value": token},
            {"Key": "MaxMinutes", "Value": str(MAX_MINUTES)},
            {"Key": "ExpiresAtUTC", "Value": deadline.isoformat()},
            {"Key": "SourceCommit", "Value": commit}]
    tag_spec = [{"ResourceType": kind, "Tags": tags}
                for kind in ("instance", "volume")]
    market = {"MarketType": "spot", "SpotOptions": {
        "SpotInstanceType": "one-time", "InstanceInterruptionBehavior": "terminate",
        "MaxPrice": MAX_PRICE}}
    mapping = [{"DeviceName": infrastructure["ami_root_device"], "Ebs": {
        "VolumeSize": VOLUME_GIB, "VolumeType": "gp3",
        "DeleteOnTermination": True, "Encrypted": True}}]
    arguments = ["ec2", "run-instances", "--region", REGION,
                 "--image-id", infrastructure["ami"],
                 "--instance-type", INSTANCE_TYPE,
                 "--subnet-id", SUBNET, "--security-group-ids", group,
                 "--iam-instance-profile", f"Name={role}",
                 "--instance-market-options", json.dumps(market),
                 "--instance-initiated-shutdown-behavior", "terminate",
                 "--metadata-options", "HttpTokens=required,HttpEndpoint=enabled",
                 "--block-device-mappings", json.dumps(mapping),
                 "--tag-specifications", json.dumps(tag_spec),
                 "--client-token", token,
                 "--user-data", f"file://{user_data}"]
    dry = command("aws", *arguments, "--dry-run", allow_error=True)
    if "DryRunOperation" not in dry.stderr:
        raise RuntimeError("exact EC2 launch dry run did not pass: " + dry.stderr)
    (output / "planned-launch.json").write_text(json.dumps({
        "client_token": token, "deadline_utc": deadline.isoformat(),
        "bucket": bucket, "security_group": group, "role": role,
        "arguments_without_user_data": arguments[:-2],
        "user_data_sha256": hashlib.sha256(
            user_data.read_bytes()).hexdigest()},
        indent=2, sort_keys=True) + "\n")
    launched = aws(*arguments)
    instance = launched["Instances"][0]
    receipt = {"schema": "n83.tier1-allocation/v1", "infrastructure": infrastructure,
               "allocation_mode": mode,
               "role": role, "bucket": bucket, "security_group": group,
               "client_token": token, "source_commit": commit,
               "instance_id": instance["InstanceId"],
               "deadline_utc": deadline.isoformat(),
               "s3_key": f"n83-tier1-online-20261016/{token}.tar.gz",
               "launch_response": launched}
    (output / "allocation.json").write_text(json.dumps(receipt, indent=2,
                                                        sort_keys=True) + "\n")
    print(json.dumps({key: receipt[key] for key in
                      ("instance_id", "client_token", "deadline_utc", "bucket",
                       "s3_key")}, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--infra-preflight", action="store_true")
    action.add_argument("--launch", action="store_true")
    action.add_argument("--launch-host-screen", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    infrastructure = current_infra()
    if args.infra_preflight:
        print(json.dumps(infrastructure, sort_keys=True))
        return 0
    if args.output is None:
        parser.error("--output is required with either launch mode")
    mode = "same_instance_selector" if args.launch_host_screen else "selected_panel"
    commit = source_preflight(mode)
    launch(args.output.resolve(), infrastructure, commit, mode)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
