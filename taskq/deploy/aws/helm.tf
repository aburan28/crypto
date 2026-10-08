locals {
  node_placement = {
    cpu = {
      nodeSelector = { "taskq/pool" = "cpu" }
      tolerations  = [{ key = "taskq/pool", operator = "Equal", value = "cpu", effect = "NoSchedule" }]
    }
    gpu = {
      nodeSelector = { "taskq/pool" = "gpu" }
      tolerations  = [{ key = "nvidia.com/gpu", operator = "Exists", effect = "NoSchedule" }]
    }
  }

  chart_pools = [for p in var.pools : merge(
    {
      name     = p.name
      queues   = p.queues
      replicas = p.replicas
      labels   = p.labels
      resources = {
        requests = merge({ cpu = p.cpu, memory = p.memory }, p.gpus > 0 ? { "nvidia.com/gpu" = p.gpus } : {})
        limits   = merge({ cpu = p.cpu, memory = p.memory }, p.gpus > 0 ? { "nvidia.com/gpu" = p.gpus } : {})
      }
    },
    local.node_placement[p.node_group],
    p.image == null ? {} : { image = p.image },
    p.keda == null ? {} : { keda = {
      minReplicas = p.keda.min_replicas
      maxReplicas = p.keda.max_replicas
      lagCount    = p.keda.lag_count
    } },
  )]

  taskq_values = {
    image = var.image
    repos = var.repos
    redis = {
      enabled           = false
      externalUrlSecret = { name = kubernetes_secret_v1.redis.metadata[0].name, key = "url" }
    }
    results = {
      create           = true
      storageClassName = kubernetes_storage_class_v1.efs.metadata[0].name
      accessMode       = "ReadWriteMany"
      storage          = var.results_storage
    }
    pools = local.chart_pools
    keda = {
      enabled = var.enable_keda
      external = {
        address        = "${aws_elasticache_replication_group.taskq.primary_endpoint_address}:6379"
        enableTLS      = true
        passwordSecret = { name = kubernetes_secret_v1.redis.metadata[0].name, key = "password" }
      }
    }
  }
}

resource "helm_release" "nvidia_device_plugin" {
  count      = var.gpu_node.enabled ? 1 : 0
  name       = "nvidia-device-plugin"
  repository = "https://nvidia.github.io/k8s-device-plugin"
  chart      = "nvidia-device-plugin"
  version    = "0.17.0"
  namespace  = "kube-system"
  values = [yamlencode({
    nodeSelector = { "taskq/pool" = "gpu" }
    tolerations  = [{ key = "nvidia.com/gpu", operator = "Exists", effect = "NoSchedule" }]
  })]
  depends_on = [module.eks]
}

resource "helm_release" "keda" {
  count            = var.enable_keda ? 1 : 0
  name             = "keda"
  repository       = "https://kedacore.github.io/charts"
  chart            = "keda"
  version          = "2.15.1"
  namespace        = "keda"
  create_namespace = true
  depends_on       = [module.eks]
}

resource "helm_release" "taskq" {
  name      = var.name
  chart     = "${path.module}/../helm/taskq"
  namespace = kubernetes_namespace_v1.taskq.metadata[0].name
  values    = [yamlencode(local.taskq_values)]
  wait      = true
  timeout   = 900

  depends_on = [helm_release.keda, helm_release.nvidia_device_plugin, aws_efs_mount_target.results]
}
