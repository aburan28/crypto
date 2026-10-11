output "cluster_name" {
  value = module.eks.cluster_name
}

output "kubeconfig_command" {
  value = "aws eks update-kubeconfig --region ${var.region} --name ${module.eks.cluster_name}"
}

output "redis_primary_endpoint" {
  value = aws_elasticache_replication_group.taskq.primary_endpoint_address
}

output "results_file_system_id" {
  value = aws_efs_file_system.results.id
}

output "submit_hint" {
  value = "kubectl -n ${var.namespace} exec deploy/${var.name}-taskq-worker-${var.pools[0].name} -- taskq stats"
}
