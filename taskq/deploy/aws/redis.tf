# ElastiCache Redis OSS: the queue. Cluster mode DISABLED, because taskq's Lua
# scripts touch several keys in one call. Multi-AZ with automatic failover,
# TLS in transit, encryption at rest, AUTH token.
#
# ElastiCache has no durable append log: a failover can drop writes that had
# not replicated. Workers mirror every result to EFS before committing it, so
# results survive; a just-submitted task may not. See README for MemoryDB.

resource "random_password" "redis_auth" {
  length  = 64
  special = false # ElastiCache forbids some specials, and the token goes in a URL
}

resource "aws_elasticache_subnet_group" "taskq" {
  name       = var.name
  subnet_ids = module.vpc.private_subnets
}

resource "aws_security_group" "redis" {
  name_prefix = "${var.name}-redis-"
  description = "taskq Redis: reachable from cluster nodes only"
  vpc_id      = module.vpc.vpc_id
  lifecycle { create_before_destroy = true }
}

resource "aws_vpc_security_group_ingress_rule" "redis_from_nodes" {
  security_group_id            = aws_security_group.redis.id
  referenced_security_group_id = module.eks.node_security_group_id
  ip_protocol                  = "tcp"
  from_port                    = 6379
  to_port                      = 6379
}

resource "aws_elasticache_replication_group" "taskq" {
  replication_group_id = var.name
  description          = "taskq queue"
  engine               = "redis"
  engine_version       = var.redis_engine_version
  node_type            = var.redis_node_type
  port                 = 6379

  num_cache_clusters         = 2
  automatic_failover_enabled = true
  multi_az_enabled           = true

  transit_encryption_enabled = true
  at_rest_encryption_enabled = true
  auth_token                 = random_password.redis_auth.result

  subnet_group_name  = aws_elasticache_subnet_group.taskq.name
  security_group_ids = [aws_security_group.redis.id]

  snapshot_retention_limit = 3
  apply_immediately        = true
}
