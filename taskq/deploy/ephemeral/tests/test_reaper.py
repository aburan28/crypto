"""Unit tests for reaper.py with a fake EC2 client (no AWS, no boto3 needed).

    python3 -m pytest taskq/deploy/ephemeral/tests/test_reaper.py
"""

import datetime as dt
import importlib.util
import pathlib
import sys
import types

# reaper.py imports boto3 at module level; a stub keeps the test dependency-free.
sys.modules.setdefault("boto3", types.SimpleNamespace(client=None))
_spec = importlib.util.spec_from_file_location(
    "reaper", pathlib.Path(__file__).resolve().parent.parent / "reaper.py"
)
reaper = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(reaper)

NOW = 1_800_000_000.0
HOUR = 3600


def instance(iid, launched_ago, deadline=None):
    tags = [{"Key": "ephemeral:managed", "Value": "true"}]
    if deadline is not None:
        tags.append({"Key": "ephemeral:deadline", "Value": str(deadline)})
    launch = dt.datetime.fromtimestamp(NOW - launched_ago, tz=dt.timezone.utc)
    return {"InstanceId": iid, "LaunchTime": launch, "Tags": tags}


class FakeEC2:
    def __init__(self, instances):
        self.instances = instances
        self.terminated = []
        self.filters = None

    def get_paginator(self, name):
        assert name == "describe_instances"
        return self

    def paginate(self, Filters):
        self.filters = Filters
        return [{"Reservations": [{"Instances": self.instances}]}]

    def terminate_instances(self, InstanceIds):
        self.terminated.extend(InstanceIds)


def run(instances):
    ec2 = FakeEC2(instances)
    result = reaper.handler(ec2=ec2, now=NOW)
    return ec2, result["terminated"]


def test_within_deadline_is_left_alone():
    _, gone = run([instance("i-ok", HOUR, deadline=NOW + HOUR)])
    assert gone == []


def test_inside_grace_is_left_alone():
    _, gone = run([instance("i-grace", 5 * HOUR, deadline=NOW - 60)])
    assert gone == []


def test_past_deadline_and_grace_is_terminated():
    _, gone = run([instance("i-late", 5 * HOUR, deadline=NOW - reaper.GRACE - 60)])
    assert gone == ["i-late"]


def test_missing_deadline_uses_max_ttl():
    young = instance("i-young", reaper.MAX_TTL - HOUR)
    old = instance("i-old", reaper.MAX_TTL + reaper.GRACE + 60)
    _, gone = run([young, old])
    assert gone == ["i-old"]


def test_garbage_deadline_uses_max_ttl():
    _, gone = run([instance("i-bad", reaper.MAX_TTL + reaper.GRACE + 60, deadline="soon")])
    assert gone == ["i-bad"]


def test_deadline_tag_cannot_exceed_max_ttl():
    far = NOW + 365 * 24 * HOUR
    _, gone = run([instance("i-greedy", reaper.MAX_TTL + reaper.GRACE + 60, deadline=far)])
    assert gone == ["i-greedy"]


def test_only_managed_instances_are_queried():
    ec2, _ = run([])
    names = {f["Name"]: f["Values"] for f in ec2.filters}
    assert names["tag:ephemeral:managed"] == ["true"]
    assert "terminated" not in names["instance-state-name"]


def test_nothing_to_do_makes_no_terminate_call():
    ec2, gone = run([instance("i-ok", HOUR, deadline=NOW + HOUR)])
    assert gone == [] and ec2.terminated == []
