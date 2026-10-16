# Standing pieces for `ephemeral run`: nothing here bills by the hour.
#
#   S3 bucket        job specs, live logs, outputs (s3://<bucket>/jobs/<id>/)
#   runner role      what a runner instance may touch: that bucket, plus SSM
#   security group   no inbound unless ssh_cidrs is set; all outbound
#   launch template  shutdown = terminate, IMDSv2, encrypted gp3 root
#   reaper Lambda    terminates any runner past its deadline tag
#   SSM parameter    /<name>/config, which the CLI reads, so any machine with
#                    the launcher policy can run jobs without this state
#   launcher policy  the rights `ephemeral` needs, to attach to users or roles

data "aws_caller_identity" "current" {}
data "aws_partition" "current" {}

locals {
  account   = data.aws_caller_identity.current.account_id
  partition = data.aws_partition.current.partition
  bucket    = "${var.name}-runs-${local.account}-${var.region}"
  managed   = "ephemeral:managed"

  # AMIs are resolved at launch time from AWS's public SSM parameters, so a
  # runner always boots the current image without a Terraform change.
  cpu_ami_param = "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64"
  arm_ami_param = "/aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-arm64"
  gpu_ami_param = "/aws/service/deeplearning/ami/x86_64/base-oss-nvidia-driver-gpu-amazon-linux-2023/latest/ami-id"

  subnet_ids = length(var.subnet_ids) > 0 ? var.subnet_ids : data.aws_subnets.default[0].ids
}

# ---- network --------------------------------------------------------------

data "aws_vpc" "default" {
  count   = length(var.subnet_ids) == 0 ? 1 : 0
  default = true
}

data "aws_subnets" "default" {
  count = length(var.subnet_ids) == 0 ? 1 : 0
  filter {
    name   = "vpc-id"
    values = [data.aws_vpc.default[0].id]
  }
  filter {
    name   = "default-for-az"
    values = ["true"]
  }
}

data "aws_subnet" "first" {
  id = local.subnet_ids[0]
}

resource "aws_security_group" "runner" {
  name        = "${var.name}-runner"
  description = "ephemeral runners: no inbound unless ssh_cidrs, all outbound"
  vpc_id      = data.aws_subnet.first.vpc_id
}

resource "aws_vpc_security_group_egress_rule" "all" {
  security_group_id = aws_security_group.runner.id
  ip_protocol       = "-1"
  cidr_ipv4         = "0.0.0.0/0"
}

resource "aws_vpc_security_group_ingress_rule" "ssh" {
  for_each          = toset(var.ssh_cidrs)
  security_group_id = aws_security_group.runner.id
  ip_protocol       = "tcp"
  from_port         = 22
  to_port           = 22
  cidr_ipv4         = each.value
}

# ---- results bucket -------------------------------------------------------

resource "aws_s3_bucket" "runs" {
  bucket = local.bucket
}

resource "aws_s3_bucket_public_access_block" "runs" {
  bucket                  = aws_s3_bucket.runs.id
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_server_side_encryption_configuration" "runs" {
  bucket = aws_s3_bucket.runs.id
  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_policy" "tls_only" {
  bucket = aws_s3_bucket.runs.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Sid       = "DenyInsecureTransport"
      Effect    = "Deny"
      Principal = "*"
      Action    = "s3:*"
      Resource  = [aws_s3_bucket.runs.arn, "${aws_s3_bucket.runs.arn}/*"]
      Condition = { Bool = { "aws:SecureTransport" = "false" } }
    }]
  })
  depends_on = [aws_s3_bucket_public_access_block.runs]
}

resource "aws_s3_bucket_lifecycle_configuration" "runs" {
  count  = var.results_expire_days > 0 ? 1 : 0
  bucket = aws_s3_bucket.runs.id
  rule {
    id     = "expire-jobs"
    status = "Enabled"
    filter {
      prefix = "jobs/"
    }
    expiration {
      days = var.results_expire_days
    }
  }
}

# ---- runner identity ------------------------------------------------------

resource "aws_iam_role" "runner" {
  name = "${var.name}-runner"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "ec2.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "runner_s3" {
  name = "runs-bucket"
  role = aws_iam_role.runner.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect   = "Allow"
        Action   = ["s3:GetObject", "s3:PutObject"]
        Resource = "${aws_s3_bucket.runs.arn}/jobs/*"
      },
      {
        Effect    = "Allow"
        Action    = "s3:ListBucket"
        Resource  = aws_s3_bucket.runs.arn
        Condition = { StringLike = { "s3:prefix" = ["jobs/*"] } }
      },
    ]
  })
}

# Session Manager: `ephemeral shell ID` without SSH keys or inbound ports.
resource "aws_iam_role_policy_attachment" "runner_ssm" {
  role       = aws_iam_role.runner.name
  policy_arn = "arn:${local.partition}:iam::aws:policy/AmazonSSMManagedInstanceCore"
}

resource "aws_iam_instance_profile" "runner" {
  name = "${var.name}-runner"
  role = aws_iam_role.runner.name
}

# ---- launch template ------------------------------------------------------

# No image here: the CLI passes resolve:ssm:<param> for the instance's
# architecture and GPU-ness. The root volume size is overridden per job.
resource "aws_launch_template" "runner" {
  name                                 = "${var.name}-runner"
  description                          = "ephemeral run: one command, then terminate"
  update_default_version               = true
  instance_initiated_shutdown_behavior = "terminate"
  key_name                             = var.key_name == "" ? null : var.key_name

  iam_instance_profile {
    arn = aws_iam_instance_profile.runner.arn
  }

  metadata_options {
    http_endpoint               = "enabled"
    http_tokens                 = "required"
    http_put_response_hop_limit = 2
    instance_metadata_tags      = "enabled"
  }

  block_device_mappings {
    device_name = "/dev/xvda"
    ebs {
      volume_size           = 100
      volume_type           = "gp3"
      encrypted             = true
      delete_on_termination = true
    }
  }

  tag_specifications {
    resource_type = "instance"
    tags          = { (local.managed) = "true" }
  }

  tag_specifications {
    resource_type = "volume"
    tags          = { (local.managed) = "true" }
  }
}

# ---- reaper ---------------------------------------------------------------

data "archive_file" "reaper" {
  type        = "zip"
  source_file = "${path.module}/reaper.py"
  output_path = "${path.module}/.build/reaper.zip"
}

resource "aws_iam_role" "reaper" {
  name = "${var.name}-reaper"
  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [{
      Effect    = "Allow"
      Principal = { Service = "lambda.amazonaws.com" }
      Action    = "sts:AssumeRole"
    }]
  })
}

resource "aws_iam_role_policy" "reaper" {
  name = "reap-runners"
  role = aws_iam_role.reaper.id
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect   = "Allow"
        Action   = "ec2:DescribeInstances"
        Resource = "*"
      },
      {
        Effect    = "Allow"
        Action    = "ec2:TerminateInstances"
        Resource  = "arn:${local.partition}:ec2:${var.region}:${local.account}:instance/*"
        Condition = { StringEquals = { "aws:ResourceTag/${local.managed}" = "true" } }
      },
      {
        Effect   = "Allow"
        Action   = ["logs:CreateLogGroup", "logs:CreateLogStream", "logs:PutLogEvents"]
        Resource = "arn:${local.partition}:logs:${var.region}:${local.account}:*"
      },
    ]
  })
}

resource "aws_lambda_function" "reaper" {
  function_name    = "${var.name}-reaper"
  role             = aws_iam_role.reaper.arn
  runtime          = "python3.12"
  handler          = "reaper.handler"
  filename         = data.archive_file.reaper.output_path
  source_code_hash = data.archive_file.reaper.output_base64sha256
  timeout          = 60
  environment {
    variables = {
      MANAGED_TAG   = local.managed
      GRACE_MINUTES = tostring(var.reaper_grace_minutes)
      MAX_TTL_HOURS = tostring(var.max_ttl_hours)
    }
  }
}

resource "aws_cloudwatch_event_rule" "reaper" {
  name                = "${var.name}-reaper"
  description         = "terminate ephemeral runners past their deadline"
  schedule_expression = var.reaper_schedule
}

resource "aws_cloudwatch_event_target" "reaper" {
  rule = aws_cloudwatch_event_rule.reaper.name
  arn  = aws_lambda_function.reaper.arn
}

resource "aws_lambda_permission" "reaper" {
  statement_id  = "AllowEventBridge"
  action        = "lambda:InvokeFunction"
  function_name = aws_lambda_function.reaper.function_name
  principal     = "events.amazonaws.com"
  source_arn    = aws_cloudwatch_event_rule.reaper.arn
}

resource "aws_cloudwatch_log_group" "reaper" {
  name              = "/aws/lambda/${aws_lambda_function.reaper.function_name}"
  retention_in_days = 30
}

# ---- config the CLI reads -------------------------------------------------

resource "aws_ssm_parameter" "config" {
  name = "/${var.name}/config"
  type = "String"
  value = jsonencode({
    region           = var.region
    bucket           = aws_s3_bucket.runs.bucket
    launch_template  = aws_launch_template.runner.id
    security_group   = aws_security_group.runner.id
    subnet_ids       = local.subnet_ids
    repo_url         = var.repo_url
    default_ttl      = var.default_ttl
    max_ttl_hours    = var.max_ttl_hours
    cpu_ami_param    = local.cpu_ami_param
    arm_ami_param    = local.arm_ami_param
    gpu_ami_param    = local.gpu_ami_param
    managed_tag      = local.managed
    instance_profile = aws_iam_instance_profile.runner.name
  })
}

# ---- launcher policy ------------------------------------------------------

resource "aws_iam_policy" "launcher" {
  name        = "${var.name}-launcher"
  description = "Rights the ephemeral CLI needs to launch, follow and kill runners"
  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Sid    = "ReadConfigAndAmis"
        Effect = "Allow"
        Action = ["ssm:GetParameter", "ssm:GetParameters"]
        Resource = [
          aws_ssm_parameter.config.arn,
          "arn:${local.partition}:ssm:${var.region}::parameter/aws/service/*",
        ]
      },
      {
        Sid    = "LaunchTaggedRunners"
        Effect = "Allow"
        Action = "ec2:RunInstances"
        Resource = [
          "arn:${local.partition}:ec2:${var.region}:${local.account}:instance/*",
          "arn:${local.partition}:ec2:${var.region}:${local.account}:volume/*",
        ]
        Condition = { StringEquals = { "aws:RequestTag/${local.managed}" = "true" } }
      },
      {
        Sid    = "LaunchFromTemplate"
        Effect = "Allow"
        Action = "ec2:RunInstances"
        Resource = [
          aws_launch_template.runner.arn,
          "arn:${local.partition}:ec2:${var.region}::image/*",
          "arn:${local.partition}:ec2:${var.region}:${local.account}:subnet/*",
          "arn:${local.partition}:ec2:${var.region}:${local.account}:network-interface/*",
          "arn:${local.partition}:ec2:${var.region}:${local.account}:security-group/*",
          "arn:${local.partition}:ec2:${var.region}:${local.account}:key-pair/*",
        ]
      },
      {
        Sid       = "TagAtLaunch"
        Effect    = "Allow"
        Action    = "ec2:CreateTags"
        Resource  = "arn:${local.partition}:ec2:${var.region}:${local.account}:*/*"
        Condition = { StringEquals = { "ec2:CreateAction" = "RunInstances" } }
      },
      {
        Sid       = "KillRunners"
        Effect    = "Allow"
        Action    = "ec2:TerminateInstances"
        Resource  = "arn:${local.partition}:ec2:${var.region}:${local.account}:instance/*"
        Condition = { StringEquals = { "aws:ResourceTag/${local.managed}" = "true" } }
      },
      {
        Sid      = "Describe"
        Effect   = "Allow"
        Action   = ["ec2:DescribeInstances", "ec2:DescribeInstanceTypes", "ec2:DescribeSpotPriceHistory", "ec2:DescribeSubnets"]
        Resource = "*"
      },
      {
        Sid      = "PassRunnerRole"
        Effect   = "Allow"
        Action   = "iam:PassRole"
        Resource = aws_iam_role.runner.arn
      },
      {
        Sid      = "RunsBucket"
        Effect   = "Allow"
        Action   = ["s3:GetObject", "s3:PutObject", "s3:ListBucket"]
        Resource = [aws_s3_bucket.runs.arn, "${aws_s3_bucket.runs.arn}/*"]
      },
      {
        Sid       = "ShellIntoRunners"
        Effect    = "Allow"
        Action    = "ssm:StartSession"
        Resource  = "arn:${local.partition}:ec2:${var.region}:${local.account}:instance/*"
        Condition = { StringEquals = { "aws:ResourceTag/${local.managed}" = "true" } }
      },
      {
        Sid    = "ShellDocument"
        Effect = "Allow"
        Action = ["ssm:StartSession", "ssm:TerminateSession"]
        Resource = [
          "arn:${local.partition}:ssm:${var.region}::document/SSM-SessionManagerRunShell",
          "arn:${local.partition}:ssm:*:*:session/$${aws:userid}-*",
        ]
      },
    ]
  })
}

resource "aws_iam_role_policy_attachment" "launcher" {
  for_each   = toset(var.launcher_role_names)
  role       = each.value
  policy_arn = aws_iam_policy.launcher.arn
}

resource "aws_iam_user_policy_attachment" "launcher" {
  for_each   = toset(var.launcher_user_names)
  user       = each.value
  policy_arn = aws_iam_policy.launcher.arn
}
