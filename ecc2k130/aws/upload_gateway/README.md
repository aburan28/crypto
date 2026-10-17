# Credential-free worker ingress

RunPod workers use a revocable campaign token instead of an AWS access key.
The canonical campaign bucket remains private and keeps S3 Block Public Access
enabled.

## Data path

```text
worker token -> HTTP API -> one-object presigned PUT -> private quarantine
                                                    -> SQS -> validator
                                                    -> canonical dp/ or ckpt/
                                                    -> existing ingester
```

The token can call only the worker API. It cannot call AWS, list the corpus,
read database credentials, delete canonical objects, or write directly to
`dp/`. The API brokers read-only access to the pinned campaign artifacts and
serializes the existing S3 slot and seed-registry compare-and-set operations.

Before an object becomes visible to the ingester, the validator checks:

- the enrolled worker, active slot owner and campaign binding;
- exact key, byte count and SHA-256, including the content-addressed key;
- v1, v2 or table-v3 record framing and a complete final record;
- the run-id and lane portions of every seed;
- no repeated seed in one upload;
- field-width bits, the canonical Frobenius-orbit representative and the
  campaign's exact `dpWeight` predicate;
- v2 witness counters against the record's iteration count;
- a gap-free, non-overlapping byte offset for each slot/stream; and
- the existing manifest fields when a strict campaign supplies one.

Rejected bodies remain only in the private quarantine bucket and expire after
two days. Accepted bodies are copied to the unchanged canonical key, so
`dp_ingest.py` and `merge.py` require no alternate corpus path.

The validator does **not** replay every v1 trail. A v1 record proves the seed is
in the leased namespace and the endpoint is a correctly encoded distinguished
orbit representative, but it does not prove that the endpoint was reached by
that seed. V2 witness checks are stronger but still do not replace the native
client's complete collision replay. Anonymous uploads are therefore not
accepted: enrollment and revocation remain part of the trust boundary.

## Deploy

The preferred serverless mode needs CloudFormation, IAM, Lambda, API Gateway,
DynamoDB, SQS and S3 permissions in `us-west-2`. It packages deterministic
Lambda code into the existing campaign bucket and creates stack
`ecc2k130-upload-gateway`:

```bash
python3 ecc2k130/aws/upload_gateway/deploy.py deploy
```

Some accounts deny Lambda service calls independently of the caller's IAM
policy. The external mode creates the same private quarantine, tables, queue
and alarm, but runs the broker and validator next to `dp_ingest.py` on one
cheap RunPod CPU pod. It leaves one narrowly scoped AWS identity on the
central service; GPU workers still receive no AWS key:

```bash
python3 ecc2k130/aws/upload_gateway/deploy.py deploy --mode external
tmpdir=$(mktemp -d)
python3 ecc2k130/aws/upload_gateway/deploy.py gateway-credentials \
  --credentials-file "$tmpdir/gateway.json"
python3 ecc2k130/aws/upload_gateway/runpod_launch.py service \
  --credentials-file "$tmpdir/gateway.json"
```

After a later gateway deployment, update an existing central pod without
putting its scoped credentials on the command line or in terminal output:

```bash
python3 ecc2k130/aws/upload_gateway/runpod_launch.py service-update POD_ID \
  --credentials-file "$tmpdir/gateway.json"
```

The launcher prints the HTTPS `apiUrl`. The credentials file is mode 0600,
is never printed, and must not be committed. Its policy is limited to the
quarantine workflow, campaign control objects, status publication, the RDS
secret, and authorizing the pod's exact egress `/32` for PostgreSQL.

The GitHub workflow can perform the same deployment with OIDC when repository
secret `ECC2K130_DEPLOY_ROLE_ARN` names the reviewed deployment role. No AWS
access key is stored in GitHub.

Enroll one worker without printing its token:

```bash
python3 ecc2k130/aws/upload_gateway/deploy.py enroll runpod-pro6000-a \
  --token-file /tmp/runpod-pro6000-a.token
```

Pass the file's value to the pod as `ECC_UPLOAD_TOKEN`, together with the
stack's `ApiUrl` as `ECC_UPLOAD_API` and the enrolled id as `ECC_WORKER_ID`.
The worker image entrypoint downloads and verifies the campaign's pinned GPU
binary, discovers the GPU through the existing supervisor, and starts the
walk. Do not pass `AWS_ACCESS_KEY_ID` or `AWS_SECRET_ACCESS_KEY`.

`runpod_launch.py` does this without putting the token on the command line or
printing it:

```bash
python3 ecc2k130/aws/upload_gateway/runpod_launch.py worker \
  --api-url https://POD_ID-8080.proxy.runpod.net \
  --worker-id runpod-pro6000-a \
  --token-file /tmp/runpod-pro6000-a.token
```

The same source bundle can be rolled into an existing broker worker without
putting its token on the command line or in terminal output:

```bash
python3 ecc2k130/aws/upload_gateway/runpod_launch.py worker-update WORKER_POD_ID \
  --api-url https://SERVICE_POD_ID-8080.proxy.runpod.net \
  --worker-id runpod-pro6000-a \
  --token-file /tmp/runpod-pro6000-a.token
```

Revoke a token immediately when a pod is retired:

```bash
python3 ecc2k130/aws/upload_gateway/deploy.py revoke \
  --token-file /tmp/runpod-pro6000-a.token
```

## Build and test

```bash
PYTHONPATH=ecc2k130/aws python3 -m unittest -v \
  ecc2k130/aws/test_upload_gateway.py \
  ecc2k130/aws/test_worker_spool.py \
  ecc2k130/aws/test_seed_registry.py
docker build -f ecc2k130/aws/upload_gateway/Dockerfile \
  -t ecc2k130-worker:local ecc2k130/aws
docker build -f ecc2k130/aws/upload_gateway/Dockerfile.gateway \
  -t ecc2k130-gateway:local ecc2k130/aws
```

On `main`, GitHub Actions publishes both tested images as
`ghcr.io/aburan28/ecc2k130-worker:<commit>` and
`ghcr.io/aburan28/ecc2k130-gateway:<commit>`, plus their `:latest` tags.
