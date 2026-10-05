# The ReadWriteMany results volume: artifacts plus the result mirror every
# worker writes before committing to Redis. Results are the one thing that
# cannot be re-run cheaply, so this is protected from `terraform destroy`.

resource "aws_efs_file_system" "results" {
  encrypted        = true
  performance_mode = "generalPurpose"
  tags             = { Name = "${var.name}-results" }
  lifecycle { prevent_destroy = true }
}

resource "aws_security_group" "efs" {
  name_prefix = "${var.name}-efs-"
  description = "taskq results EFS: NFS from cluster nodes only"
  vpc_id      = module.vpc.vpc_id
  lifecycle { create_before_destroy = true }
}

resource "aws_vpc_security_group_ingress_rule" "efs_from_nodes" {
  security_group_id            = aws_security_group.efs.id
  referenced_security_group_id = module.eks.node_security_group_id
  ip_protocol                  = "tcp"
  from_port                    = 2049
  to_port                      = 2049
}

resource "aws_efs_mount_target" "results" {
  # Keyed by AZ, which is known at plan time; subnet IDs are not until the
  # VPC exists, and for_each refuses unknown keys.
  for_each        = { for i, az in local.azs : az => i }
  file_system_id  = aws_efs_file_system.results.id
  subnet_id       = module.vpc.private_subnets[each.value]
  security_groups = [aws_security_group.efs.id]
}
