# Offline plan with a mocked AWS provider: checks the bootstrap plans and that
# the state bucket and lock table carry the protections the main stack relies on.
mock_provider "aws" {
  mock_data "aws_caller_identity" {
    defaults = { account_id = "123456789012" }
  }
  mock_data "aws_iam_policy_document" {
    defaults = { json = "{\"Version\":\"2012-10-17\",\"Statement\":[]}" }
  }
}
mock_provider "random" {}

run "bootstrap_plans" {
  command = plan
  assert {
    condition     = aws_dynamodb_table.locks.hash_key == "LockID" && aws_dynamodb_table.locks.name == "taskq-tflock"
    error_message = "the S3 backend locks on a LockID hash key"
  }
  assert {
    condition     = aws_s3_bucket_versioning.state.versioning_configuration[0].status == "Enabled"
    error_message = "state bucket must be versioned"
  }
  assert {
    condition = alltrue([
      aws_s3_bucket_public_access_block.state.block_public_acls,
      aws_s3_bucket_public_access_block.state.block_public_policy,
      aws_s3_bucket_public_access_block.state.ignore_public_acls,
      aws_s3_bucket_public_access_block.state.restrict_public_buckets,
    ])
    error_message = "state bucket must block all public access"
  }
  assert {
    condition     = one(one(aws_s3_bucket_server_side_encryption_configuration.state.rule).apply_server_side_encryption_by_default).sse_algorithm == "aws:kms"
    error_message = "state bucket must be KMS-encrypted"
  }
}
