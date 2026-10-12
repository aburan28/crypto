locals {
  # Static CPU manager: a Guaranteed pod with integer CPUs owns its cores, which
  # is what makes two measurements comparable. It requires reserved CPU > 0.
  static_cpu_manager = [{
    content_type = "application/node.eks.aws"
    content      = <<-EOT
      ---
      apiVersion: node.eks.aws/v1alpha1
      kind: NodeConfig
      spec:
        kubelet:
          config:
            cpuManagerPolicy: static
            kubeReserved: {cpu: 500m, memory: 1Gi}
            systemReserved: {cpu: 300m, memory: 512Mi}
    EOT
  }]

  worker_node_groups = merge(
    {
      cpu = {
        ami_type              = "AL2023_x86_64_STANDARD"
        instance_types        = [var.cpu_node.instance_type]
        min_size              = var.cpu_node.min_size
        desired_size          = var.cpu_node.desired_size
        max_size              = var.cpu_node.max_size
        labels                = { "taskq/pool" = "cpu" }
        taints                = { pool = { key = "taskq/pool", value = "cpu", effect = "NO_SCHEDULE" } }
        cloudinit_pre_nodeadm = local.static_cpu_manager
      }
    },
    var.gpu_node.enabled ? {
      gpu = {
        ami_type              = "AL2023_x86_64_NVIDIA" # drivers in the AMI; device plugin below
        instance_types        = [var.gpu_node.instance_type]
        min_size              = var.gpu_node.min_size
        desired_size          = var.gpu_node.desired_size
        max_size              = var.gpu_node.max_size
        labels                = { "taskq/pool" = "gpu" }
        taints                = { gpu = { key = "nvidia.com/gpu", value = "present", effect = "NO_SCHEDULE" } }
        cloudinit_pre_nodeadm = local.static_cpu_manager
      }
    } : {}
  )
}

module "eks" {
  source  = "terraform-aws-modules/eks/aws"
  version = "~> 20.24"

  cluster_name    = var.name
  cluster_version = var.cluster_version

  cluster_endpoint_public_access           = true
  cluster_endpoint_public_access_cidrs     = var.cluster_endpoint_public_access_cidrs
  enable_cluster_creator_admin_permissions = true

  vpc_id     = module.vpc.vpc_id
  subnet_ids = module.vpc.private_subnets

  cluster_addons = {
    coredns    = {}
    kube-proxy = {}
    vpc-cni    = {}
    aws-efs-csi-driver = {
      service_account_role_arn = module.efs_csi_irsa.iam_role_arn
    }
  }

  eks_managed_node_groups = merge(
    {
      # KEDA, the CSI controller, the device plugin and the MCP endpoint.
      system = {
        ami_type       = "AL2023_x86_64_STANDARD"
        instance_types = ["t3.large"]
        min_size       = 1
        desired_size   = 2
        max_size       = 3
      }
    },
    local.worker_node_groups
  )
}

module "efs_csi_irsa" {
  source  = "terraform-aws-modules/iam/aws//modules/iam-role-for-service-accounts-eks"
  version = "~> 5.44"

  role_name_prefix      = "${var.name}-efs-csi-"
  attach_efs_csi_policy = true
  oidc_providers = {
    main = {
      provider_arn               = module.eks.oidc_provider_arn
      namespace_service_accounts = ["kube-system:efs-csi-controller-sa"]
    }
  }
}
