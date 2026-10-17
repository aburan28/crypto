terraform {
  required_version = ">= 1.6"
  required_providers {
    aws     = { source = "hashicorp/aws", version = "~> 5.70" }
    archive = { source = "hashicorp/archive", version = "~> 2.6" }
  }

  # Nothing secret lives in this state, so the default local backend is fine
  # for a single operator. To share it, add `backend "s3" {}` and point it at
  # the bucket that ../aws/bootstrap creates.
}

provider "aws" {
  region = var.region
  default_tags {
    tags = merge({ Project = "ephemeral-runner", managed-by = "terraform" }, var.tags)
  }
}
