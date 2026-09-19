#!/usr/bin/env python3
"""Tests for audit_userdata.py, with a fake EC2 in place of the account."""

import base64
import unittest

import audit_userdata as audit

REAL_PREAMBLE = """#!/bin/bash
install -d -m 700 /root/.aws /var/lib/ecc2k130
cat > /root/.aws/credentials <<'AWSCREDS'
[default]
aws_access_key_id=AKIAEXAMPLEEXAMPLE12
aws_secret_access_key=notarealsecretnotarealsecretnotarealsec
AWSCREDS
"""

SESSION_PREAMBLE = """#!/bin/bash
cat > /root/.aws/credentials <<'CREDS'
[default]
aws_access_key_id=ASIAEXAMPLEEXAMPLE34
aws_secret_access_key=alsonotarealsecret
aws_session_token=IQoJb3JpZ2luX2VjEAAaCXVzLXdlc3QtMiJ
CREDS
"""

CLEAN = """#!/bin/bash
export AWS_DEFAULT_REGION=us-west-2
aws s3 cp "s3://bucket/campaign.json" campaign.json
"""


def b64(text):
    return base64.b64encode(text.encode()).decode()


class FakePaginator:
    def __init__(self, pages):
        self.pages = pages

    def paginate(self, **kwargs):
        for page in self.pages:
            yield page


class FakeEc2:
    """Only the four calls the audit makes, shaped like botocore's."""

    def __init__(self, templates=(), versions=None, instances=(), user_data=None):
        self.templates = list(templates)
        self.versions = versions or {}
        self.instances = list(instances)
        self.user_data = user_data or {}
        self.attribute_calls = []

    def get_paginator(self, name):
        if name == "describe_launch_templates":
            return FakePaginator([{"LaunchTemplates":
                                   [{"LaunchTemplateName": t} for t in self.templates]}])
        if name == "describe_launch_template_versions":
            outer = self

            class ByName:
                def paginate(self, LaunchTemplateName=None, **kwargs):
                    yield {"LaunchTemplateVersions": outer.versions.get(LaunchTemplateName, [])}
            return ByName()
        if name == "describe_instances":
            return FakePaginator([{"Reservations": [{"Instances": self.instances}]}])
        raise AssertionError("unexpected paginator %s" % name)

    def describe_instance_attribute(self, InstanceId=None, Attribute=None):
        self.attribute_calls.append((InstanceId, Attribute))
        value = self.user_data.get(InstanceId)
        return {"UserData": {"Value": value} if value else {}}


class Findings(unittest.TestCase):
    def test_it_finds_a_long_lived_key_and_its_secret(self):
        self.assertEqual(audit.findings(REAL_PREAMBLE),
                         ["access key id", "secret key assignment"])

    def test_it_finds_temporary_session_credentials(self):
        self.assertEqual(audit.findings(SESSION_PREAMBLE),
                         ["access key id", "secret key assignment", "session token assignment"])

    def test_a_worker_bootstrap_without_credentials_is_clean(self):
        self.assertEqual(audit.findings(CLEAN), [])

    def test_the_words_alone_are_not_a_finding(self):
        """A comment about credentials is not a credential."""
        self.assertEqual(audit.findings(
            "# the instance profile supplies aws credentials through IMDS\n"
            "# do not write an aws_secret_access_key here\n"), [])

    def test_it_does_not_report_the_secret_itself(self):
        rows = audit.audit_region(
            FakeEc2(templates=["t"],
                    versions={"t": [{"VersionNumber": 1,
                                     "LaunchTemplateData": {"UserData": b64(REAL_PREAMBLE)}}]}),
            "us-west-2")
        blob = repr(rows)
        self.assertNotIn("notarealsecret", blob)
        self.assertNotIn("AKIAEXAMPLEEXAMPLE12", blob)


class AuditRegion(unittest.TestCase):
    def test_it_names_the_template_version_and_that_it_lacks_a_profile(self):
        ec2 = FakeEc2(
            templates=["ecc2k130-worker"],
            versions={"ecc2k130-worker": [
                {"VersionNumber": 40,
                 "LaunchTemplateData": {"UserData": b64(REAL_PREAMBLE)}},
                {"VersionNumber": 42,
                 "LaunchTemplateData": {"UserData": b64(CLEAN),
                                        "IamInstanceProfile": {"Name": "ecc2k130-worker"}}},
            ]})
        rows = audit.audit_region(ec2, "us-west-2")
        self.assertEqual(len(rows), 1, "only the tainted version is reported")
        self.assertEqual(rows[0]["id"], "ecc2k130-worker version 40")
        self.assertIsNone(rows[0]["profile"])

    def test_it_checks_stopped_instances_too(self):
        """The old ingest host was stopped, and still readable."""
        ec2 = FakeEc2(
            instances=[{"InstanceId": "i-clean"}, {"InstanceId": "i-tainted"}],
            user_data={"i-clean": b64(CLEAN), "i-tainted": b64(SESSION_PREAMBLE)})
        rows = audit.audit_region(ec2, "us-east-2")
        self.assertEqual([r["id"] for r in rows], ["i-tainted"])
        self.assertEqual(
            sorted(ec2.attribute_calls),
            [("i-clean", "userData"), ("i-tainted", "userData")])

    def test_a_clean_account_reports_nothing(self):
        ec2 = FakeEc2(
            templates=["ecc2k130-worker"],
            versions={"ecc2k130-worker": [
                {"VersionNumber": 42,
                 "LaunchTemplateData": {"UserData": b64(CLEAN),
                                        "IamInstanceProfile": {"Name": "ecc2k130-worker"}}}]},
            instances=[{"InstanceId": "i-1"}], user_data={"i-1": b64(CLEAN)})
        self.assertEqual(audit.audit_region(ec2, "us-west-2"), [])

    def test_an_instance_with_no_user_data_is_not_an_error(self):
        ec2 = FakeEc2(instances=[{"InstanceId": "i-none"}])
        self.assertEqual(audit.audit_region(ec2, "us-west-2"), [])

    def test_user_data_that_is_not_base64_is_still_read(self):
        ec2 = FakeEc2(instances=[{"InstanceId": "i-plain"}],
                      user_data={"i-plain": REAL_PREAMBLE})
        self.assertEqual([r["id"] for r in audit.audit_region(ec2, "us-west-2")], ["i-plain"])


if __name__ == "__main__":
    unittest.main()
