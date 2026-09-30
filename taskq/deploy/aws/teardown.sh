#!/usr/bin/env bash
# Remove what provision.sh created. The EFS results volume is KEPT unless
# DELETE_RESULTS=1: results are the one thing that cannot be re-run cheaply.
set -euo pipefail
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
CLUSTER=${CLUSTER:-taskq}
CACHE_ID=${CACHE_ID:-taskq}
NAMESPACE=${NAMESPACE:-taskq}
RELEASE=${RELEASE:-rq}

helm uninstall "$RELEASE" -n "$NAMESPACE" 2>/dev/null || true
if aws elasticache delete-replication-group --replication-group-id "$CACHE_ID" --no-retain-primary-cluster 2>/dev/null; then
    aws elasticache wait replication-group-deleted --replication-group-id "$CACHE_ID"
fi
aws elasticache delete-cache-subnet-group --cache-subnet-group-name "$CACHE_ID" 2>/dev/null || true
if [ "${DELETE_RESULTS:-0}" = 1 ]; then
    FS=$(aws efs describe-file-systems --query "FileSystems[?Name=='$CLUSTER-results'].FileSystemId | [0]" --output text)
    if [ -n "$FS" ] && [ "$FS" != None ]; then
        for mt in $(aws efs describe-mount-targets --file-system-id "$FS" --query 'MountTargets[].MountTargetId' --output text); do
            aws efs delete-mount-target --mount-target-id "$mt"
        done
        while [ "$(aws efs describe-mount-targets --file-system-id "$FS" --query 'length(MountTargets)')" != 0 ]; do sleep 5; done
        aws efs delete-file-system --file-system-id "$FS"
    fi
else
    echo "keeping EFS $CLUSTER-results (DELETE_RESULTS=1 to delete; eksctl delete will then fail on its SG until you do)"
fi
eksctl delete cluster --name "$CLUSTER" --wait
