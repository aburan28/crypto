#!/usr/bin/env bash
set -euo pipefail

REGION=${AWS_REGION:-us-west-2}
NAME=${NAME:-crypto-g7e-dev}
INSTANCE_TYPE=${INSTANCE_TYPE:-g7e.2xlarge}
KEY_NAME=${KEY_NAME:-meow34}
ROOT_GB=${ROOT_GB:-500}
SECURITY_GROUP_NAME=${SECURITY_GROUP_NAME:-crypto-g7e-dev-ssh}
BOOTSTRAP=${BOOTSTRAP:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/bootstrap.sh}
SSH_USER=${SSH_USER:-ubuntu}

aws_ec2() { aws --region "$REGION" ec2 "$@"; }

instance_id() {
  local id
  if ! id=$(aws_ec2 describe-instances \
      --filters "Name=tag:Name,Values=$NAME" "Name=instance-state-name,Values=pending,running,stopping,stopped" \
      --query 'Reservations[].Instances[] | sort_by(@,&LaunchTime)[-1].InstanceId' \
      --output text); then
    return 1
  fi
  if [ "$id" != "None" ]; then
    printf '%s\n' "$id"
  fi
}

public_ip() {
  local id=${1:-$(instance_id)}
  [ -n "$id" ] || return 1
  aws_ec2 describe-instances --instance-ids "$id" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' --output text
}

ensure_sg() {
  local vpc_id sg_id ssh_cidr
  sg_id=${1:-}
  if [ -z "$sg_id" ]; then
    vpc_id=$(aws_ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    [ "$vpc_id" != "None" ] || { echo "No default VPC; set SG_ID to an existing security group." >&2; exit 2; }

    sg_id=$(aws_ec2 describe-security-groups --filters "Name=group-name,Values=$SECURITY_GROUP_NAME" "Name=vpc-id,Values=$vpc_id" --query 'SecurityGroups[0].GroupId' --output text 2>/dev/null || true)
    if [ -z "$sg_id" ] || [ "$sg_id" = "None" ]; then
      sg_id=$(aws_ec2 create-security-group --group-name "$SECURITY_GROUP_NAME" --description "SSH for durable G7e developer host" --vpc-id "$vpc_id" --query GroupId --output text)
    fi
  fi

  ssh_cidr=${SSH_CIDR:-}
  if [ -z "$ssh_cidr" ]; then
    ssh_cidr="$(curl -fsS https://checkip.amazonaws.com | tr -d '\n')/32"
  fi
  aws_ec2 authorize-security-group-ingress --group-id "$sg_id" --protocol tcp --port 22 --cidr "$ssh_cidr" >/dev/null 2>&1 || true
  printf '%s\n' "$sg_id"
}

latest_dlami() {
  aws_ec2 describe-images --owners amazon \
    --filters 'Name=name,Values=Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 24.04)*' 'Name=state,Values=available' 'Name=architecture,Values=x86_64' \
    --query 'reverse(sort_by(Images,&CreationDate))[0].ImageId' --output text
}

cmd_up() {
  local id state ami root_dev sg_id instance_sg_id
  id=$(instance_id)
  if [ -n "$id" ]; then
    state=$(aws_ec2 describe-instances --instance-ids "$id" --query 'Reservations[0].Instances[0].State.Name' --output text)
    instance_sg_id=$(aws_ec2 describe-instances --instance-ids "$id" \
      --query 'Reservations[0].Instances[0].SecurityGroups[0].GroupId' --output text)
    ensure_sg "$instance_sg_id" >/dev/null
    case "$state" in
      stopped)
        echo "Starting existing durable host $id ..."
        aws_ec2 start-instances --instance-ids "$id" >/dev/null
        ;;
      stopping)
        echo "Waiting for existing durable host $id to stop ..."
        aws_ec2 wait instance-stopped --instance-ids "$id"
        echo "Starting existing durable host $id ..."
        aws_ec2 start-instances --instance-ids "$id" >/dev/null
        ;;
      *)
        echo "Existing host: $id ($state)"
        ;;
    esac
    aws_ec2 wait instance-running --instance-ids "$id"
    echo "Instance: $id"
    echo "Public IP: $(public_ip "$id")"
    return
  fi

  [ -f "$BOOTSTRAP" ] || { echo "Missing bootstrap: $BOOTSTRAP" >&2; exit 2; }
  ami=${AMI_ID:-$(latest_dlami)}
  [ "$ami" != "None" ] || { echo "Could not resolve a G7e-compatible Ubuntu DLAMI; set AMI_ID." >&2; exit 2; }
  root_dev=$(aws_ec2 describe-images --image-ids "$ami" --query 'Images[0].RootDeviceName' --output text)
  sg_id=${SG_ID:-$(ensure_sg)}

  echo "Launching $INSTANCE_TYPE on-demand in $REGION (AMI $ami, key $KEY_NAME) ..."
  id=$(aws_ec2 run-instances \
    --image-id "$ami" \
    --instance-type "$INSTANCE_TYPE" \
    --key-name "$KEY_NAME" \
    --security-group-ids "$sg_id" \
    --user-data "file://$BOOTSTRAP" \
    --instance-initiated-shutdown-behavior stop \
    --block-device-mappings "DeviceName=$root_dev,Ebs={VolumeSize=$ROOT_GB,VolumeType=gp3,Encrypted=true,DeleteOnTermination=false}" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$NAME},{Key=Purpose,Value=durable-gpu-devhost},{Key=Lifecycle,Value=on-demand},{Key=CostGuardManaged,Value=true},{Key=CostGuardMonthlyBudget,Value=${MONTHLY_BUDGET_USD:-5000}}]" \
    --metadata-options 'HttpTokens=required,HttpEndpoint=enabled' \
    --query 'Instances[0].InstanceId' --output text)

  # Protect against accidental termination. Stop/start is the normal lifecycle.
  aws_ec2 modify-instance-attribute --instance-id "$id" --disable-api-termination >/dev/null
  aws_ec2 wait instance-running --instance-ids "$id"
  echo "Instance: $id"
  echo "Public IP: $(public_ip "$id")"
  echo "SSH: ssh -i meow34.pem $SSH_USER@$(public_ip "$id")"
}

cmd_stop() {
  local id
  id=$(instance_id)
  [ -n "$id" ] || { echo "No host found"; return; }
  aws_ec2 stop-instances --instance-ids "$id" >/dev/null
  echo "Stopping $id; EBS and installed tools persist."
}

cmd_status() {
  local id
  id=$(instance_id)
  [ -n "$id" ] || { echo "No host found"; return; }
  aws_ec2 describe-instances --instance-ids "$id" \
    --query 'Reservations[0].Instances[0].{InstanceId:InstanceId,State:State.Name,Type:InstanceType,AZ:Placement.AvailabilityZone,PublicIp:PublicIpAddress,LaunchTime:LaunchTime}' \
    --output table
}

cmd_ssh() {
  local ip
  ip=$(public_ip)
  [ -n "$ip" ] && [ "$ip" != "None" ] || { echo "Instance has no public IP / is not running." >&2; exit 2; }
  chmod 600 meow34.pem
  exec ssh -i meow34.pem "$SSH_USER@$ip"
}

cmd_reinstall() {
  local ip
  ip=$(public_ip)
  [ -n "$ip" ] && [ "$ip" != "None" ] || { echo "Instance must be running." >&2; exit 2; }
  chmod 600 meow34.pem
  ssh -i meow34.pem "$SSH_USER@$ip" 'sudo /usr/local/bin/refresh-agent-tools'
}

case ${1:-status} in
  up|start) cmd_up ;;
  stop|down) cmd_stop ;;
  status) cmd_status ;;
  ip) public_ip ;;
  ssh) cmd_ssh ;;
  reinstall|refresh) cmd_reinstall ;;
  *) echo "usage: $0 {up|stop|status|ip|ssh|reinstall}" >&2; exit 2 ;;
esac
