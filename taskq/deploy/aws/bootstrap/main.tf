# One-time bootstrap: the S3 bucket and DynamoDB lock table that hold the
# main stack's (../) Terraform state. This module keeps its own state LOCALLY
# (terraform.tfstate here, gitignored): it holds no secrets, and it cannot
# store its state in the bucket it is about to create.
#
#   terraform init && terraform apply
#   terraform output -raw backend_hcl > ../backend.hcl
#
# Both resources are protected from `terraform destroy`: losing them loses the
# record of what the main stack created.

terraform {
  required_version = ">= 1.6"
  required_providers {
    aws    = { source = "hashicorp/aws", version = "~> 5.70" }
    random = { source = "hashicorp/random", version = "~> 3.6" }
  }
}

variable "region" {
  type    = string
  default = "us-west-2"
}

variable "name" {
  description = "Prefix; must match the main stack's var.name for the default state key."
  type        = string
  default     = "taskq"
}

provider "aws" {
  region = var.region
  default_tags {
    tags = { project = "taskq", managed-by = "terraform", component = "state-bootstrap" }
  }
}

data "aws_caller_identity" "current" {}

# Bucket names are global: account id plus a random suffix avoids collisions.
resource "random_id" "suffix" {
  byte_length = 3
}

resource "aws_s3_bucket" "state" {
  bucket = "${var.name}-tfstate-${data.aws_caller_identity.current.account_id}-${random_id.suffix.hex}"
  lifecycle { prevent_destroy = true }
}

resource "aws_s3_bucket_versioning" "state" {
  bucket = aws_s3_bucket.state.id
  versioning_configuration { status = "Enabled" } # recover a corrupted or overwritten state
}

resource "aws_s3_bucket_server_side_encryption_configuration" "state" {
  bucket = aws_s3_bucket.state.id
  rule {
    apply_server_side_encryption_by_default { sse_algorithm = "aws:kms" }
    bucket_key_enabled = true
  }
}

resource "aws_s3_bucket_public_access_block" "state" {
  bucket                  = aws_s3_bucket.state.id
  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

resource "aws_s3_bucket_ownership_controls" "state" {
  bucket = aws_s3_bucket.state.id
  rule { object_ownership = "BucketOwnerEnforced" }
}

data "aws_iam_policy_document" "tls_only" {
  statement {
    sid       = "DenyInsecureTransport"
    effect    = "Deny"
    actions   = ["s3:*"]
    resources = [aws_s3_bucket.state.arn, "${aws_s3_bucket.state.arn}/*"]
    principals {
      type        = "*"
      identifiers = ["*"]
    }
    condition {
      test     = "Bool"
      variable = "aws:SecureTransport"
      values   = ["false"]
    }
  }
}

resource "aws_s3_bucket_policy" "state" {
  bucket     = aws_s3_bucket.state.id
  policy     = data.aws_iam_policy_document.tls_only.json
  depends_on = [aws_s3_bucket_public_access_block.state]
}

resource "aws_s3_bucket_lifecycle_configuration" "state" {
  bucket = aws_s3_bucket.state.id
  rule {
    id     = "expire-old-state-versions"
    status = "Enabled"
    filter {}
    noncurrent_version_expiration { noncurrent_days = 90 }
  }
  depends_on = [aws_s3_bucket_versioning.state]
}

resource "aws_dynamodb_table" "locks" {
  name         = "${var.name}-tflock"
  billing_mode = "PAY_PER_REQUEST"
  hash_key     = "LockID" # the key name the S3 backend uses
  attribute {
    name = "LockID"
    type = "S"
  }
  server_side_encryption { enabled = true }
  point_in_time_recovery { enabled = true }
  lifecycle { prevent_destroy = true }
}

output "bucket" {
  value = aws_s3_bucket.state.id
}

output "lock_table" {
  value = aws_dynamodb_table.locks.name
}

output "backend_hcl" {
  description = "terraform output -raw backend_hcl > ../backend.hcl"
  value       = <<-EOT
    bucket         = "${aws_s3_bucket.state.id}"
    key            = "${var.name}/${var.region}/terraform.tfstate"
    region         = "${var.region}"
    encrypt        = true
    dynamodb_table = "${aws_dynamodb_table.locks.name}"
  EOT
}
