variable "region" {
  description = "AWS region. us-west-2 matches ecc2k130/aws."
  type        = string
  default     = "us-west-2"
}

variable "name" {
  description = "Prefix for every resource (cluster, cache, file system)."
  type        = string
  default     = "taskq"
}

variable "tags" {
  description = "Extra tags on every resource."
  type        = map(string)
  default     = {}
}

variable "vpc_cidr" {
  type    = string
  default = "10.42.0.0/16"
}

variable "cluster_version" {
  type    = string
  default = "1.31"
}

variable "cluster_endpoint_public_access_cidrs" {
  description = "Who may reach the Kubernetes API. Terraform itself needs it; narrow it to your own addresses."
  type        = list(string)
  default     = ["0.0.0.0/0"]
}

# ---- node groups ------------------------------------------------------------
# One worker pod per node: size nodes for one pool's resources plus reserved.
# EKS managed node groups ignore desired_size after creation (the module keeps
# it out of the diff so autoscalers can own it); change a running group's size
# with `aws eks update-nodegroup-config` or install an autoscaler.

variable "cpu_node" {
  type = object({
    instance_type = string
    min_size      = number
    desired_size  = number
    max_size      = number
  })
  default = { instance_type = "c7i.2xlarge", min_size = 0, desired_size = 2, max_size = 20 }
}

variable "gpu_node" {
  description = "g7e = RTX PRO 6000 Blackwell (sm_120). desired_size 0: a GPU costs money while idle."
  type = object({
    enabled       = bool
    instance_type = string
    min_size      = number
    desired_size  = number
    max_size      = number
  })
  default = { enabled = true, instance_type = "g7e.2xlarge", min_size = 0, desired_size = 0, max_size = 4 }
}

# ---- Redis ------------------------------------------------------------------

variable "redis_node_type" {
  type    = string
  default = "cache.r7g.large"
}

variable "redis_engine_version" {
  description = "Redis OSS 7.x: taskq needs XAUTOCLAIM (>= 6.2) and Lua."
  type        = string
  default     = "7.1"
}

# ---- workloads --------------------------------------------------------------

variable "image" {
  description = "Default worker image. Pin a sha tag so every result's TASKQ_IMAGE names what ran."
  type        = object({ repository = string, tag = string })
  default     = { repository = "ghcr.io/aburan28/taskq-worker", tag = "latest" }
}

variable "repos" {
  description = "Worker repo allowlist: alias -> clone URL."
  type        = map(string)
  default = {
    crypto                = "https://github.com/aburan28/crypto.git"
    crypto-autoresearcher = "https://github.com/aburan28/crypto-autoresearcher.git"
    cryptanalysis         = "https://github.com/aburan28/cryptanalysis.git"
  }
}

variable "pools" {
  description = <<-EOT
    Worker pools, one Deployment each. node_group is "cpu" or "gpu" and decides
    the node selector and tolerations. cpu/memory are both request and limit
    (Guaranteed QoS, so the static CPU manager gives the pod exclusive cores).
    image overrides var.image (e.g. the cryptanalysis toolchain images).
  EOT
  type = list(object({
    name       = string
    queues     = list(string)
    replicas   = number
    node_group = string
    cpu        = string
    memory     = string
    gpus       = optional(number, 0)
    labels     = optional(map(string), {})
    image      = optional(object({ repository = string, tag = string }))
    # KEDA scaler for this pool's queues; only used when enable_keda = true.
    keda = optional(object({
      min_replicas = optional(number, 0)
      max_replicas = optional(number, 10)
      lag_count    = optional(number, 1)
    }))
  }))
  default = [
    { name = "cpu", queues = ["cpu"], replicas = 2, node_group = "cpu",
    cpu = "6", memory = "12Gi", labels = { pool = "cpu-c7i" } },
  ]
  validation {
    condition     = alltrue([for p in var.pools : contains(["cpu", "gpu"], p.node_group)])
    error_message = "pools[*].node_group must be \"cpu\" or \"gpu\"."
  }
}

variable "namespace" {
  type    = string
  default = "taskq"
}

variable "enable_keda" {
  description = "KEDA redis-streams scalers. Pods only; without a node autoscaler extra pods sit Pending."
  type        = bool
  default     = false
}

variable "results_storage" {
  description = "Nominal PVC size; EFS is elastic and bills for what is stored."
  type        = string
  default     = "100Gi"
}
