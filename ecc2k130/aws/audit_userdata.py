#!/usr/bin/env python3
"""Find AWS credentials embedded in EC2 user-data.

User-data is not a secret store. Every process on the instance can read it
through IMDS, so can anything that gets a shell in a container on it, and so can
any principal in the account holding `ec2:DescribeLaunchTemplateVersions` or
`ec2:DescribeInstanceAttribute` -- which is a much larger set than the people
who are supposed to hold the campaign's keys.

This campaign put a long-lived key there for weeks. `infra.sh` prefers an
instance profile and only falls back to `WORKER_AWS_ACCESS_KEY_ID` /
`WORKER_AWS_SECRET_ACCESS_KEY` when the profile cannot be created, and for a
while the caller could not create one, so every worker launched carrying the
static credentials of an IAM user with AdministratorAccess. Six other launchers
here still have the same fallback. A fallback nobody audits is a fallback that
becomes the default, so this makes the audit one command:

    python3 audit_userdata.py                     # this account, usual regions
    python3 audit_userdata.py --regions us-west-2 # one region
    python3 audit_userdata.py --quiet             # exit status only

Exit status is 0 when nothing is found and 1 when something is, so it can gate a
deployment. It reports where the credential is, never what it is.
"""

import argparse
import base64
import re
import sys

DEFAULT_REGIONS = ("us-west-2", "us-east-1", "us-east-2")

# Long-lived keys (AKIA), temporary session keys (ASIA), and the assignment
# forms that carry the secret half. The key id alone is not secret; it is
# reported because it is the reliable marker that the secret is beside it.
PATTERNS = (
    ("access key id", re.compile(r"\b(?:AKIA|ASIA)[0-9A-Z]{12,}\b")),
    ("secret key assignment", re.compile(r"(?i)\baws_secret_access_key\b\s*[=:]")),
    ("session token assignment", re.compile(r"(?i)\baws_session_token\b\s*[=:]")),
)


def findings(text):
    """What kinds of credential material appear in this user-data."""
    return sorted({name for name, rx in PATTERNS if rx.search(text)})


def decode(value):
    if not value:
        return ""
    try:
        return base64.b64decode(value).decode("utf-8", "replace")
    except Exception:
        # Console-entered user-data is not always base64; read it as text.
        return value if isinstance(value, str) else ""


def audit_region(ec2, region):
    """Every launch-template version and every instance that still has one."""
    rows = []
    templates = []
    paginator = ec2.get_paginator("describe_launch_templates")
    for page in paginator.paginate():
        templates.extend(t["LaunchTemplateName"] for t in page["LaunchTemplates"])
    for name in templates:
        versions = []
        vp = ec2.get_paginator("describe_launch_template_versions")
        for page in vp.paginate(LaunchTemplateName=name):
            versions.extend(page["LaunchTemplateVersions"])
        for version in versions:
            data = version.get("LaunchTemplateData", {})
            hits = findings(decode(data.get("UserData")))
            if hits:
                rows.append({
                    "region": region,
                    "kind": "launch-template",
                    "id": "%s version %s" % (name, version["VersionNumber"]),
                    "profile": (data.get("IamInstanceProfile") or {}).get("Name"),
                    "found": hits,
                })
    ip = ec2.get_paginator("describe_instances")
    instances = []
    for page in ip.paginate(Filters=[{"Name": "instance-state-name",
                                      "Values": ["pending", "running", "stopping", "stopped"]}]):
        for reservation in page["Reservations"]:
            instances.extend(reservation["Instances"])
    for instance in instances:
        iid = instance["InstanceId"]
        attribute = ec2.describe_instance_attribute(InstanceId=iid, Attribute="userData")
        hits = findings(decode((attribute.get("UserData") or {}).get("Value")))
        if hits:
            rows.append({
                "region": region,
                "kind": "instance",
                "id": iid,
                "profile": (instance.get("IamInstanceProfile") or {}).get("Arn"),
                "found": hits,
            })
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", nargs="+", default=list(DEFAULT_REGIONS))
    parser.add_argument("--quiet", action="store_true", help="exit status only")
    args = parser.parse_args(argv)

    import boto3

    rows = []
    for region in args.regions:
        rows.extend(audit_region(boto3.client("ec2", region_name=region), region))

    if not args.quiet:
        for row in rows:
            print("%s  %-15s %-40s %s%s" % (
                row["region"], row["kind"], row["id"], ", ".join(row["found"]),
                "" if row["profile"] else "  (and no instance profile)"))
        print("%d place(s) with credentials in user-data across %s"
              % (len(rows), ", ".join(args.regions)))
    return 1 if rows else 0


if __name__ == "__main__":
    sys.exit(main())
