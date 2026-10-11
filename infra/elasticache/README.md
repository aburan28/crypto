# Gröbner-basis ElastiCache

This stack creates a private, two-node Redis OSS replication group for
index-calculus workers. It enables TLS, at-rest encryption, automatic failover,
and a Secrets Manager-generated auth token. The cache is scoped to the worker
security group; it is not publicly reachable.

The Rust solver remains fully local when `IC_GROEBNER_CACHE_URL` is unset.
When it is set, matrix-F4 results are cached by a versioned BLAKE3 key over the
complete Boolean system, variable count, and reduction degree. Redis failures
are treated as cache misses, so the cache cannot change a solver result.

## Provision

Use at least two private subnets in different Availability Zones and pass the
security group attached to the worker instances:

```bash
aws cloudformation deploy \
  --region "$AWS_REGION" \
  --stack-name ic-groebner-cache \
  --template-file infra/elasticache/groebner-cache.yaml \
  --parameter-overrides \
    VpcId=vpc-0123456789abcdef0 \
    SubnetIds='subnet-0123456789abcdef0,subnet-fedcba9876543210' \
    ApplicationSecurityGroupId=sg-0123456789abcdef0
```

Replace every example ID. The worker security group must be in the same VPC
as the cache subnets. The command creates the auth secret; do not put that
secret in a `.env` file committed to the repository.

Retrieve the endpoint and secret only on the worker host:

```bash
export IC_CACHE_ENDPOINT="$(
  aws cloudformation describe-stacks \
    --stack-name ic-groebner-cache \
    --query 'Stacks[0].Outputs[?OutputKey==`PrimaryEndpointAddress`].OutputValue' \
    --output text
)"
export IC_CACHE_PORT="$(
  aws cloudformation describe-stacks \
    --stack-name ic-groebner-cache \
    --query 'Stacks[0].Outputs[?OutputKey==`PrimaryEndpointPort`].OutputValue' \
    --output text
)"
export IC_CACHE_SECRET_ARN="$(
  aws cloudformation describe-stacks \
    --stack-name ic-groebner-cache \
    --query 'Stacks[0].Outputs[?OutputKey==`AuthSecretArn`].OutputValue' \
    --output text
)"
export IC_CACHE_AUTH_TOKEN="$(
  aws secretsmanager get-secret-value \
    --secret-id "$IC_CACHE_SECRET_ARN" \
    --query SecretString \
    --output text |
    python3 -c 'import json, sys; print(json.load(sys.stdin)["password"])'
)"
export IC_GROEBNER_CACHE_URL="rediss://:${IC_CACHE_AUTH_TOKEN}@${IC_CACHE_ENDPOINT}:${IC_CACHE_PORT}/0"
export IC_GROEBNER_CACHE_TTL_SECS=604800
```

The `rediss://` scheme is required because the stack requires in-transit
encryption. The Redis client validates the server certificate using the host
system's native TLS trust store.

## Tune or disable

Optional settings:

```bash
export IC_GROEBNER_CACHE_NAMESPACE=crypto:ic:groebner:v1
export IC_GROEBNER_CACHE_MAX_BYTES=$((16 * 1024 * 1024))
```

Unset `IC_GROEBNER_CACHE_URL` to disable remote caching. Cached values are
recomputable artifacts, not source-of-truth research data; the stack retains a
final snapshot on deletion for operational recovery, but entries can be
evicted at any time.

## Remove

```bash
aws cloudformation delete-stack \
  --region "$AWS_REGION" \
  --stack-name ic-groebner-cache
```

The cache and its generated secret are intentionally infrastructure-managed.
Review the retained snapshot and secret before deleting them separately.
