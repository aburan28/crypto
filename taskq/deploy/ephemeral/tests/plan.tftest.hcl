# Offline applies against mocked providers: no AWS account needed. They check that
# the configuration applies and that the safety properties are wired in, not
# that AWS would accept it.
mock_provider "aws" {
  mock_data "aws_caller_identity" {
    defaults = { account_id = "123456789012" }
  }
  mock_data "aws_partition" {
    defaults = { partition = "aws" }
  }
  mock_data "aws_subnets" {
    defaults = { ids = ["subnet-a", "subnet-b", "subnet-c"] }
  }
  mock_data "aws_subnet" {
    defaults = { vpc_id = "vpc-0test" }
  }
  mock_resource "aws_iam_role" {
    defaults = { arn = "arn:aws:iam::123456789012:role/mock" }
  }
  mock_resource "aws_iam_instance_profile" {
    defaults = { arn = "arn:aws:iam::123456789012:instance-profile/mock" }
  }
  mock_resource "aws_s3_bucket" {
    defaults = { arn = "arn:aws:s3:::mock" }
  }
  mock_resource "aws_lambda_function" {
    defaults = { arn = "arn:aws:lambda:us-east-1:123456789012:function:mock" }
  }
  mock_resource "aws_ssm_parameter" {
    defaults = { arn = "arn:aws:ssm:us-east-1:123456789012:parameter/ephemeral/config" }
  }
  mock_resource "aws_cloudwatch_event_rule" {
    defaults = { arn = "arn:aws:events:us-east-1:123456789012:rule/ephemeral-reaper" }
  }
  mock_resource "aws_iam_policy" {
    defaults = { arn = "arn:aws:iam::123456789012:policy/ephemeral-launcher" }
  }
  mock_resource "aws_launch_template" {
    defaults = { arn = "arn:aws:ec2:us-east-1:123456789012:launch-template/lt-mock", id = "lt-mock" }
  }
}

run "defaults" {
  command = apply

  assert {
    condition     = aws_launch_template.runner.instance_initiated_shutdown_behavior == "terminate"
    error_message = "a runner's shutdown must terminate it"
  }
  assert {
    condition     = aws_launch_template.runner.metadata_options[0].http_tokens == "required"
    error_message = "runners must require IMDSv2"
  }
  assert {
    condition     = aws_launch_template.runner.metadata_options[0].instance_metadata_tags == "enabled"
    error_message = "runners read their tags from IMDS"
  }
  assert {
    condition     = aws_launch_template.runner.block_device_mappings[0].ebs[0].encrypted == "true"
    error_message = "root volumes must be encrypted"
  }
  assert {
    condition     = length(aws_vpc_security_group_ingress_rule.ssh) == 0
    error_message = "no inbound rule unless ssh_cidrs is set"
  }
  assert {
    condition     = aws_s3_bucket.runs.bucket == "ephemeral-runs-123456789012-us-east-1"
    error_message = "bucket name is <name>-runs-<account>-<region>"
  }
  assert {
    condition     = jsondecode(aws_ssm_parameter.config.value).subnet_ids == ["subnet-a", "subnet-b", "subnet-c"]
    error_message = "config lists the default VPC's subnets"
  }
  assert {
    condition     = aws_cloudwatch_event_rule.reaper.schedule_expression == "rate(10 minutes)"
    error_message = "the reaper runs on a schedule"
  }
  assert {
    condition     = length(aws_s3_bucket_lifecycle_configuration.runs) == 0
    error_message = "results are kept by default"
  }
}

run "custom_network_and_ssh" {
  command = apply
  variables {
    subnet_ids          = ["subnet-x"]
    ssh_cidrs           = ["203.0.113.4/32"]
    key_name            = ""
    results_expire_days = 30
    launcher_role_names = ["InstanceRole"]
  }

  assert {
    condition     = jsondecode(aws_ssm_parameter.config.value).subnet_ids == ["subnet-x"]
    error_message = "explicit subnets replace the default VPC"
  }
  assert {
    condition     = length(aws_vpc_security_group_ingress_rule.ssh) == 1
    error_message = "one SSH rule per CIDR"
  }
  assert {
    condition     = aws_launch_template.runner.key_name == null
    error_message = "an empty key_name launches without a key"
  }
  assert {
    condition     = length(aws_s3_bucket_lifecycle_configuration.runs) == 1
    error_message = "results_expire_days adds a lifecycle rule"
  }
  assert {
    condition     = length(aws_iam_role_policy_attachment.launcher) == 1
    error_message = "launcher policy is attached to the named roles"
  }
}
