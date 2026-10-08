# Offline plan with mocked providers: no AWS account, no cluster. It checks
# that the configuration plans (types, for_each keys known at plan time,
# references between modules), not that AWS would accept it.
mock_provider "aws" {
  mock_data "aws_availability_zones" {
    defaults = { names = ["us-west-2a", "us-west-2b", "us-west-2c"] }
  }
  mock_data "aws_iam_policy_document" {
    defaults = { json = "{\"Version\":\"2012-10-17\",\"Statement\":[]}" }
  }
  mock_data "aws_partition" {
    defaults = { partition = "aws", dns_suffix = "amazonaws.com" }
  }
  mock_data "aws_caller_identity" {
    defaults = { account_id = "123456789012", arn = "arn:aws:iam::123456789012:user/test" }
  }
  mock_data "aws_iam_session_context" {
    defaults = { issuer_arn = "arn:aws:iam::123456789012:user/test" }
  }
  mock_resource "aws_eks_cluster" {
    defaults = {
      endpoint              = "https://example.eks.amazonaws.com"
      certificate_authority = [{ data = "Y2E=" }]
      identity              = [{ oidc = [{ issuer = "https://oidc.eks.us-west-2.amazonaws.com/id/EXAMPLE" }] }]
    }
  }
  mock_resource "aws_iam_role" {
    defaults = { arn = "arn:aws:iam::123456789012:role/mock" }
  }
  mock_resource "aws_iam_policy" {
    defaults = { arn = "arn:aws:iam::123456789012:policy/mock" }
  }
  mock_resource "aws_iam_openid_connect_provider" {
    defaults = { arn = "arn:aws:iam::123456789012:oidc-provider/oidc.eks.us-west-2.amazonaws.com/id/EXAMPLE" }
  }
}
mock_provider "kubernetes" {}
mock_provider "helm" {}
mock_provider "random" {}
mock_provider "tls" {}
mock_provider "time" {}
mock_provider "cloudinit" {}
mock_provider "null" {}

variables {
  cluster_endpoint_public_access_cidrs = ["203.0.113.0/24"]
  enable_keda                          = true
  pools = [
    { name = "cpu", queues = ["cpu"], replicas = 2, node_group = "cpu", cpu = "6", memory = "12Gi",
    keda = { max_replicas = 20 } },
    { name = "ca-gpu-sm120", queues = ["ca-gpu-sm120"], replicas = 0, node_group = "gpu",
      cpu  = "6", memory = "48Gi", gpus = 1,
    image = { repository = "ghcr.io/aburan28/cryptanalysis-taskq-gpu", tag = "latest" } },
  ]
}

run "plans_with_gpu_and_keda" {
  command = plan
  assert {
    condition     = length(aws_efs_mount_target.results) == 3
    error_message = "one EFS mount target per AZ"
  }
  assert {
    condition     = length(helm_release.nvidia_device_plugin) == 1 && length(helm_release.keda) == 1
    error_message = "device plugin and KEDA must be installed when enabled"
  }
  assert {
    condition     = aws_elasticache_replication_group.taskq.transit_encryption_enabled && aws_elasticache_replication_group.taskq.multi_az_enabled
    error_message = "Redis must be TLS and Multi-AZ"
  }
  assert {
    condition     = local.chart_pools[0].keda.maxReplicas == 20 && local.chart_pools[1].resources.limits["nvidia.com/gpu"] == 1
    error_message = "pool keda block and GPU limits must reach the chart"
  }
}

run "cpu_only_without_keda" {
  command = plan
  variables {
    enable_keda = false
    gpu_node    = { enabled = false, instance_type = "g7e.2xlarge", min_size = 0, desired_size = 0, max_size = 4 }
    pools       = [{ name = "cpu", queues = ["cpu"], replicas = 1, node_group = "cpu", cpu = "6", memory = "12Gi" }]
  }
  assert {
    condition     = length(helm_release.nvidia_device_plugin) == 0 && length(helm_release.keda) == 0
    error_message = "no GPU or KEDA components when disabled"
  }
}

run "rejects_unknown_node_group" {
  command = plan
  variables {
    pools = [{ name = "x", queues = ["x"], replicas = 1, node_group = "tpu", cpu = "1", memory = "1Gi" }]
  }
  expect_failures = [var.pools]
}
