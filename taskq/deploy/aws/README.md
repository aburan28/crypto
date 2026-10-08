# taskq on AWS: Terraform for EKS + ElastiCache + EFS

This directory is one Terraform root module. It creates the whole stack and
installs taskq onto it:

| file | creates |
|---|---|
| `network.tf` | VPC across 3 AZs, private and public subnets, one NAT gateway |
| `eks.tf` | EKS 1.31 with three managed node groups: `system` (t3.large), `cpu` (c7i.2xlarge) and `gpu` (g7e.2xlarge, RTX PRO 6000 Blackwell, sm_120). The worker groups run the kubelet's **static CPU manager**, so a Guaranteed worker pod owns its cores and its timings are comparable. Also the EFS CSI driver with its IAM role |
| `redis.tf` | ElastiCache Redis 7.1: cluster mode disabled, a primary and a replica, Multi-AZ failover, TLS, encryption at rest, a generated AUTH token, reachable only from the cluster nodes' security group |
| `efs.tf` | encrypted EFS results volume with a mount target per AZ, reachable only from the nodes; `prevent_destroy` |
| `kubernetes.tf` | the `taskq` namespace, the `taskq-redis` Secret (URL and password) and the `efs-sc` StorageClass |
| `helm.tf` | the `deploy/helm/taskq` chart wired to all of the above, the NVIDIA device plugin when the GPU group is enabled, and KEDA when `enable_keda` is set |

## Use

```sh
cd taskq/deploy/aws
(cd bootstrap && terraform init && terraform apply \
   && terraform output -raw backend_hcl > ../backend.hcl)   # one-time: state bucket + lock table
cp terraform.tfvars.example terraform.tfvars  # pools, sizes, API CIDRs
terraform init -backend-config=backend.hcl
terraform plan -out tfplan
terraform apply tfplan
$(terraform output -raw kubeconfig_command)
kubectl -n taskq get pods
```

Prerequisites:

* AWS credentials with rights to create a VPC, EKS, EC2, IAM roles, ElastiCache and EFS.
* The `aws` CLI on the machine that runs Terraform. The Kubernetes and Helm providers get their cluster token from `aws eks get-token`.
* Terraform 1.6 or later.

**`bootstrap/` creates the state backend once.** It makes a versioned,
KMS-encrypted S3 bucket and a DynamoDB lock table. The bucket blocks all public
access and refuses non-TLS requests. Both resources are protected from
`terraform destroy`. The bootstrap keeps its own small state locally, in
`bootstrap/terraform.tfstate` (gitignored, no secrets). Keep that file, or
`terraform import` the two resources later. To use a bucket you already have
instead, fill in `backend.hcl.example` by hand.

**State holds the Redis AUTH token.** Keep state in the encrypted S3 backend
and never commit it. `.gitignore` here excludes state, `terraform.tfvars` and
`backend.hcl`.

## Pools

`pools` in `terraform.tfvars` defines one worker Deployment per entry:

* `node_group` is `cpu` or `gpu`. It picks the node selector and tolerations.
* `cpu`/`memory` are both the request and the limit, which is Guaranteed QoS.
* `gpus` adds `nvidia.com/gpu`.
* `image` overrides the default worker image, for example with the cryptanalysis toolchain images.
* `keda` sets that pool's scaler bounds when `enable_keda = true`.

The example tfvars has three pools:

* `cpu`, for crypto and crypto-autoresearcher jobs;
* `ca-cpu`, for cryptanalysis `ca` jobs;
* `ca-gpu-sm120`, which starts at 0 replicas.

Size nodes for **one worker pod each**. A c7i.2xlarge has 8 vCPU and 16 GiB.
After kube-reserved and system-reserved, a 6 CPU / 12 GiB pod fills it.

## Scaling and cost

* **No node autoscaler is installed.** The EKS module also ignores
  `desired_size` after a node group is created, so changing tfvars does not
  resize a live group. Resize it directly, then match the pool replicas:

  ```sh
  aws eks update-nodegroup-config --cluster-name taskq \
      --nodegroup-name <cpu group name from the console or `aws eks list-nodegroups`> \
      --scaling-config minSize=0,maxSize=20,desiredSize=6
  ```

* **KEDA scales pods, not nodes.** Without Karpenter or Cluster Autoscaler,
  pods KEDA adds sit Pending until nodes exist. It is off by default.
* **The GPU group starts at 0 nodes.** To use it:
  1. Scale the group to 1.
  2. Set `ca-gpu-sm120` replicas to 1 and apply.
  3. Scale back to 0 when the queue is empty: a g7e bills every hour it runs.
* **Costs while the stack is up:**
  * the EKS control plane;
  * the NAT gateway;
  * 2 t3.large system nodes, plus the c7i nodes;
  * 2 ElastiCache nodes;
  * EFS storage.

## Durability

ElastiCache has no durable append log. A primary failure followed by a
failover can lose writes that had not replicated: a just-submitted task, or a
just-committed state change. Results are not lost, because workers write every
result to EFS (`/results/results/<task>/attempt-N-fence-F.json`) *before*
committing it to Redis.

If the queue itself must survive failover, use **Amazon MemoryDB**. It speaks
the Redis protocol and has a durable Multi-AZ transaction log. To use it:

* Replace `aws_elasticache_replication_group` with `aws_memorydb_cluster`, with one shard.
* Set the chart's `namespace: "{taskq}"` so every key hashes to one slot.

This module does not do that switch yet.

## Teardown

`terraform destroy` removes everything except the EFS results volume, which
has `prevent_destroy`. To delete the results as well, remove that line first.

## Checks

```sh
terraform fmt -check -recursive
terraform init -backend=false && terraform validate
terraform test     # offline plans with mocked providers: no AWS account needed
(cd bootstrap && terraform init -input=false && terraform validate && terraform test)
```

`tests/plan.tftest.hcl` plans the whole stack against mocked providers in
three cases: with GPU and KEDA, CPU-only, and an invalid pool, which must be
refused. It checks that the configuration plans: types, `for_each` keys known
at plan time, and wiring between modules. It does not check that AWS would
accept every value. CI runs all three checks.

## Images

The default image is `ghcr.io/aburan28/taskq-worker`. The cryptanalysis
images are `ghcr.io/aburan28/cryptanalysis-taskq-{cpu,gpu}`. All three are
public: they inherit visibility from their public repositories and pull
anonymously, so no pull secret is needed. Pin `tag` to a commit sha in
production, so every result's `TASKQ_IMAGE` names what really ran.
