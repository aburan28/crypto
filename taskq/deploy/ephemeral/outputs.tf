output "bucket" {
  value = aws_s3_bucket.runs.bucket
}

output "config_parameter" {
  description = "Set EPHEMERAL_CONFIG to this if var.name is not the default."
  value       = aws_ssm_parameter.config.name
}

output "launcher_policy_arn" {
  description = "Attach to any user or role that should run `ephemeral run`."
  value       = aws_iam_policy.launcher.arn
}

output "runner_role_arn" {
  value = aws_iam_role.runner.arn
}

output "launch_template_id" {
  value = aws_launch_template.runner.id
}
