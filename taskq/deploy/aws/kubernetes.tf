resource "kubernetes_namespace_v1" "taskq" {
  metadata { name = var.namespace }
  depends_on = [module.eks]
}

# The only copy of the AUTH token outside Terraform state. Workers read the
# URL; KEDA reads the bare password.
resource "kubernetes_secret_v1" "redis" {
  metadata {
    name      = "taskq-redis"
    namespace = kubernetes_namespace_v1.taskq.metadata[0].name
  }
  data = {
    url      = "rediss://:${random_password.redis_auth.result}@${aws_elasticache_replication_group.taskq.primary_endpoint_address}:6379/0"
    password = random_password.redis_auth.result
  }
}

resource "kubernetes_storage_class_v1" "efs" {
  metadata { name = "efs-sc" }
  storage_provisioner = "efs.csi.aws.com"
  parameters = {
    provisioningMode = "efs-ap"
    fileSystemId     = aws_efs_file_system.results.id
    directoryPerms   = "700"
    uid              = "10001" # the worker image's user
    gid              = "10001"
  }
  depends_on = [module.eks, aws_efs_mount_target.results]
}
