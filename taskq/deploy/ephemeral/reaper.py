"""Terminate ephemeral runners that outlived their deadline.

The runner arms its own shutdown timer as the first line of user data, so this
only fires when that never happened: the instance hung before cloud-init, the
timer was cancelled by hand, or the OS wedged. A runner's deadline is its
`ephemeral:deadline` tag (Unix seconds), set by the CLI at launch. A runner
with no readable deadline gets launch time + MAX_TTL_HOURS.
"""

import os
import time

import boto3

MANAGED_TAG = os.environ.get("MANAGED_TAG", "ephemeral:managed")
DEADLINE_TAG = "ephemeral:deadline"
GRACE = int(os.environ.get("GRACE_MINUTES", "15")) * 60
MAX_TTL = int(float(os.environ.get("MAX_TTL_HOURS", "48")) * 3600)


def overdue(instance, now):
    """Return the number of seconds past deadline + grace, or None."""
    tags = {t["Key"]: t["Value"] for t in instance.get("Tags", [])}
    launched = instance["LaunchTime"].timestamp()
    try:
        deadline = float(tags[DEADLINE_TAG])
    except (KeyError, ValueError):
        deadline = launched + MAX_TTL
    # Never let a tag extend a runner beyond the hard ceiling.
    deadline = min(deadline, launched + MAX_TTL)
    late = now - (deadline + GRACE)
    return late if late > 0 else None


def handler(event=None, context=None, ec2=None, now=None):
    ec2 = ec2 or boto3.client("ec2")
    now = now if now is not None else time.time()
    doomed = []
    pages = ec2.get_paginator("describe_instances").paginate(
        Filters=[
            {"Name": f"tag:{MANAGED_TAG}", "Values": ["true"]},
            {"Name": "instance-state-name", "Values": ["pending", "running", "stopping", "stopped"]},
        ]
    )
    for page in pages:
        for reservation in page["Reservations"]:
            for instance in reservation["Instances"]:
                late = overdue(instance, now)
                if late is not None:
                    doomed.append(instance["InstanceId"])
                    print(f"reaping {instance['InstanceId']}: {late / 60:.0f} min past deadline + grace")
    if doomed:
        ec2.terminate_instances(InstanceIds=doomed)
    return {"terminated": doomed}
