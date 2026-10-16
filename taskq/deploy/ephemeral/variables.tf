variable "region" {
  description = "Region the runners launch in. us-east-1 is where the ECC2K-130 GPU quota lives."
  type        = string
  default     = "us-east-1"
}

variable "name" {
  description = "Prefix for every resource, and the SSM config parameter /<name>/config."
  type        = string
  default     = "ephemeral"
}

variable "repo_url" {
  description = "Git URL the runners clone. Public HTTPS needs no credentials."
  type        = string
  default     = "https://github.com/aburan28/crypto.git"
}

variable "subnet_ids" {
  description = "Subnets to launch into. Empty means the default VPC's default subnets (one per AZ). They must give instances a route to the internet (public IP or NAT)."
  type        = list(string)
  default     = []
}

variable "key_name" {
  description = "EC2 key pair for SSH (AGENTS.md section 9: meow34). Empty launches without a key; `ephemeral shell` uses SSM instead."
  type        = string
  default     = "meow34"
}

variable "ssh_cidrs" {
  description = "CIDRs allowed to SSH in. Empty opens no inbound port at all."
  type        = list(string)
  default     = []
}

variable "default_ttl" {
  description = "Hard lifetime of a runner when `ephemeral run` is given no --ttl (systemd time span)."
  type        = string
  default     = "4h"
}

variable "max_ttl_hours" {
  description = "The CLI refuses a --ttl longer than this, and the reaper uses it for runners whose deadline tag is missing or unreadable."
  type        = number
  default     = 48
}

variable "reaper_grace_minutes" {
  description = "The reaper terminates a runner this long after its deadline tag, if the on-box shutdown timer has not already."
  type        = number
  default     = 15
}

variable "reaper_schedule" {
  description = "EventBridge schedule for the reaper Lambda."
  type        = string
  default     = "rate(10 minutes)"
}

variable "results_expire_days" {
  description = "Delete job logs and outputs from S3 after this many days. 0 keeps them forever."
  type        = number
  default     = 0
}

variable "launcher_role_names" {
  description = "Existing IAM roles to attach the launcher policy to (for example the dev box's InstanceRole)."
  type        = list(string)
  default     = []
}

variable "launcher_user_names" {
  description = "Existing IAM users to attach the launcher policy to."
  type        = list(string)
  default     = []
}

variable "tags" {
  description = "Extra tags for every resource."
  type        = map(string)
  default     = {}
}
