#!/usr/bin/env bash
# Stand up taskq on AWS: EKS (eksctl), ElastiCache Redis, EFS results volume,
# and the Helm release. Idempotent: re-running it skips what exists.
#
#   AWS_DEFAULT_REGION  us-west-2     CLUSTER   taskq       NAMESPACE  taskq
#   CACHE_ID            taskq         CACHE_NODE_TYPE  cache.r7g.large
#   RELEASE             rq            INSTALL_KEDA     0
#
# Needs: aws (v2), eksctl, kubectl, helm, openssl. Costs money: see README.md.
set -euo pipefail
cd "$(dirname "$0")"

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
CLUSTER=${CLUSTER:-taskq}
NAMESPACE=${NAMESPACE:-taskq}
CACHE_ID=${CACHE_ID:-taskq}
CACHE_NODE_TYPE=${CACHE_NODE_TYPE:-cache.r7g.large}
RELEASE=${RELEASE:-rq}
INSTALL_KEDA=${INSTALL_KEDA:-0}

for tool in aws eksctl kubectl helm openssl; do
    command -v "$tool" >/dev/null || { echo "missing: $tool" >&2; exit 1; }
done
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
echo "account $ACCOUNT, region $AWS_DEFAULT_REGION, cluster $CLUSTER"

# ---- EKS --------------------------------------------------------------------
if ! aws eks describe-cluster --name "$CLUSTER" >/dev/null 2>&1; then
    sed "s/^  name: taskq$/  name: $CLUSTER/; s/^  region: us-west-2$/  region: $AWS_DEFAULT_REGION/" \
        cluster.yaml > /tmp/taskq-cluster.yaml
    eksctl create cluster -f /tmp/taskq-cluster.yaml
else
    echo "cluster $CLUSTER exists"
    aws eks update-kubeconfig --name "$CLUSTER" >/dev/null
fi
VPC=$(aws eks describe-cluster --name "$CLUSTER" --query cluster.resourcesVpcConfig.vpcId --output text)
CIDR=$(aws ec2 describe-vpcs --vpc-ids "$VPC" --query 'Vpcs[0].CidrBlock' --output text)
# eksctl tags its private subnets with the internal-elb role.
SUBNETS=$(aws ec2 describe-subnets --filters "Name=vpc-id,Values=$VPC" \
    "Name=tag:kubernetes.io/role/internal-elb,Values=1" --query 'Subnets[].SubnetId' --output text)
[ -n "$SUBNETS" ] || { echo "no private subnets found in $VPC" >&2; exit 1; }

# Security group: VPC-internal only, Redis and NFS. No public ingress anywhere.
sg() {
    local name=$1 port=$2 id
    id=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$name" "Name=vpc-id,Values=$VPC" \
         --query 'SecurityGroups[0].GroupId' --output text 2>/dev/null || true)
    if [ -z "$id" ] || [ "$id" = None ]; then
        id=$(aws ec2 create-security-group --group-name "$name" --vpc-id "$VPC" \
             --description "taskq $name (VPC-internal)" --query GroupId --output text)
        aws ec2 authorize-security-group-ingress --group-id "$id" --protocol tcp \
            --port "$port" --cidr "$CIDR" >/dev/null
    fi
    echo "$id"
}

# ---- ElastiCache ------------------------------------------------------------
# Cluster mode DISABLED: taskq's Lua scripts touch several keys at once.
# Multi-AZ with a replica and automatic failover, TLS in transit, AUTH token.
# ElastiCache has no durable log: a failover can drop the last writes. Workers
# mirror every result to EFS before committing, so no result is lost; for a
# durable queue too, use MemoryDB (see README.md).
REDIS_SG=$(sg "$CLUSTER-redis" 6379)
if ! aws elasticache describe-cache-subnet-groups --cache-subnet-group-name "$CACHE_ID" >/dev/null 2>&1; then
    # shellcheck disable=SC2086
    aws elasticache create-cache-subnet-group --cache-subnet-group-name "$CACHE_ID" \
        --cache-subnet-group-description "taskq" --subnet-ids $SUBNETS >/dev/null
fi
if ! aws elasticache describe-replication-groups --replication-group-id "$CACHE_ID" >/dev/null 2>&1; then
    TOKEN=$(openssl rand -hex 32)
    aws elasticache create-replication-group --replication-group-id "$CACHE_ID" \
        --replication-group-description "taskq queue" --engine redis --engine-version 7.1 \
        --cache-node-type "$CACHE_NODE_TYPE" --num-cache-clusters 2 \
        --automatic-failover-enabled --multi-az-enabled \
        --transit-encryption-enabled --at-rest-encryption-enabled --auth-token "$TOKEN" \
        --cache-subnet-group-name "$CACHE_ID" --security-group-ids "$REDIS_SG" >/dev/null
    echo "creating ElastiCache $CACHE_ID (about 10 minutes)"
    aws elasticache wait replication-group-available --replication-group-id "$CACHE_ID"
    NEW_TOKEN=1
else
    echo "ElastiCache $CACHE_ID exists"
    NEW_TOKEN=0
fi
REDIS_HOST=$(aws elasticache describe-replication-groups --replication-group-id "$CACHE_ID" \
    --query 'ReplicationGroups[0].NodeGroups[0].PrimaryEndpoint.Address' --output text)

# ---- EFS (ReadWriteMany results) --------------------------------------------
EFS_SG=$(sg "$CLUSTER-efs" 2049)
FS=$(aws efs describe-file-systems --query "FileSystems[?Name=='$CLUSTER-results'].FileSystemId | [0]" --output text)
if [ -z "$FS" ] || [ "$FS" = None ]; then
    FS=$(aws efs create-file-system --encrypted --performance-mode generalPurpose \
         --tags "Key=Name,Value=$CLUSTER-results" --query FileSystemId --output text)
    until [ "$(aws efs describe-file-systems --file-system-id "$FS" --query 'FileSystems[0].LifeCycleState' --output text)" = available ]; do sleep 5; done
    for s in $SUBNETS; do
        aws efs create-mount-target --file-system-id "$FS" --subnet-id "$s" --security-groups "$EFS_SG" >/dev/null
    done
    echo "created EFS $FS"
fi

# ---- Kubernetes objects -----------------------------------------------------
kubectl get namespace "$NAMESPACE" >/dev/null 2>&1 || kubectl create namespace "$NAMESPACE"
kubectl apply -f - <<YAML
apiVersion: storage.k8s.io/v1
kind: StorageClass
metadata: {name: efs-sc}
provisioner: efs.csi.aws.com
parameters:
  provisioningMode: efs-ap
  fileSystemId: $FS
  directoryPerms: "700"
  uid: "10001"
  gid: "10001"
YAML
if [ "$NEW_TOKEN" = 1 ]; then
    # The token exists only here and in this Secret; it is never printed.
    kubectl -n "$NAMESPACE" create secret generic taskq-redis \
        --from-literal=url="rediss://:$TOKEN@$REDIS_HOST:6379/0" \
        --from-literal=password="$TOKEN" --dry-run=client -o yaml | kubectl apply -f -
elif ! kubectl -n "$NAMESPACE" get secret taskq-redis >/dev/null 2>&1; then
    echo "ElastiCache exists but secret $NAMESPACE/taskq-redis does not; recreate it with the AUTH token" >&2
    exit 1
fi

KEDA_ARGS=()
if [ "$INSTALL_KEDA" = 1 ]; then
    helm repo add kedacore https://kedacore.github.io/charts >/dev/null
    helm upgrade --install keda kedacore/keda -n keda --create-namespace --wait
    KEDA_ARGS=(--set keda.enabled=true)
fi
helm upgrade --install "$RELEASE" ../helm/taskq -n "$NAMESPACE" -f values-aws.yaml \
    --set "keda.external.address=$REDIS_HOST:6379" "${KEDA_ARGS[@]}" --wait --timeout 15m

kubectl -n "$NAMESPACE" get pods -o wide
echo
echo "Submit from a pod in the cluster, or port-forward through one:"
echo "  kubectl -n $NAMESPACE exec deploy/$RELEASE-taskq-worker-cpu -- taskq stats"
