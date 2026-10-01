# taskq on AWS: EKS + ElastiCache

```sh
cd taskq/deploy/aws
./provision.sh               # EKS, ElastiCache, EFS, secrets, Helm release (idempotent)
INSTALL_KEDA=1 ./provision.sh   # optionally add KEDA scalers
./teardown.sh                # everything except the EFS results volume
```

`provision.sh` needs valid AWS credentials, plus `aws` (v2), `eksctl`,
`kubectl`, `helm` and `openssl`. The region defaults to `us-west-2`, the same
as `ecc2k130/aws`. Every name can be overridden by environment variable; see
the script header.

| component | what | why |
|---|---|---|
| EKS `taskq` | `cluster.yaml`: a `system` managed group, plus self-managed `cpu-c7i` (c7i.2xlarge) and `gpu-g7e` (g7e.2xlarge, RTX PRO 6000 Blackwell) groups | Self-managed groups run the kubelet's **static CPU manager**, so a Guaranteed worker pod owns its cores and its timings are comparable |
| ElastiCache `taskq` | Redis 7.1, cluster mode **disabled**, primary + replica, Multi-AZ failover, TLS, AUTH token, `cache.r7g.large` | The queue. Cluster mode is off because taskq's Lua scripts touch several keys at once |
| EFS `taskq-results` | encrypted, with a mount target per private subnet, as StorageClass `efs-sc` (uid 10001) | The ReadWriteMany artifacts and result-mirror volume that every worker writes |
| Helm release `rq` | `values-aws.yaml` over `deploy/helm/taskq` | CPU pool of 2 pods, each 6 CPU / 12 GiB, one per c7i node. The GPU pool is ready but commented out |

Nothing gets a public endpoint. Redis and EFS accept connections only from
inside the VPC. To reach the queue, use `kubectl exec` into a worker, or
port-forward through a pod. The AUTH token is generated here, written only to
the Secret `taskq/taskq-redis`, and never printed.

## Durability: read this before you trust the queue

ElastiCache has no durable append log. A primary failure followed by a
failover can lose the last writes that had not replicated yet. For taskq
that means a just-submitted task, or a just-committed state change, can be
lost. It does **not** mean a result can be lost, because workers write every
result to EFS (`/results/results/<task>/attempt-N-fence-F.json`) *before*
committing it to Redis. If the queue itself must survive failover, swap in
**Amazon MemoryDB**. It speaks the Redis protocol with a durable Multi-AZ
transaction log. A single-shard MemoryDB cluster works with a plain client.
Set `namespace: "{taskq}"` in the chart values so every key hashes to one
slot. `provision.sh` does not create MemoryDB yet.

## Scaling and cost

No node autoscaler is installed, so KEDA is off by default. Extra pods would
sit Pending. Scale nodes and pods together:

```sh
eksctl scale nodegroup --cluster taskq --name cpu-c7i --nodes 6
helm upgrade rq ../helm/taskq -n taskq -f values-aws.yaml --reuse-values --set pools[0].replicas=6
```

The GPU group starts at 0 nodes. To use it, publish the cryptanalysis GPU
image, uncomment the `gpu-sm120` pool, scale `gpu-g7e` to 1 and upgrade.
Scale it back to 0 when the queue is empty, because a g7e costs money every
hour it runs. Other costs while the stack is up: the EKS control plane, two
t3.large system nodes, the c7i nodes, the two ElastiCache nodes, and EFS
storage.

## Images

The chart pulls `ghcr.io/aburan28/taskq-worker`. The `taskq image` workflow
publishes it on each merge. GHCR packages start private, so either make the
package public or add an `imagePullSecrets` entry. The cryptanalysis images
come from that repo's `cloud/taskq/`.
