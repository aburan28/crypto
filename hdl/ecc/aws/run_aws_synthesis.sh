#!/usr/bin/env bash
#
# run_aws_synthesis.sh -- synthesise hdl/ecc on an AWS build instance and
# bring the utilisation and timing reports back.
#
# WHY THIS EXISTS.  docs/ecc_fpga_cost_model.md estimates the datapath at
# ~256 DSPs and ~300 MHz, and those estimates carry the whole FPGA-vs-GPU
# comparison.  This script replaces them with numbers from Vivado.
#
# WHAT IT DOES NOT DO.  It does not build an Amazon FPGA Image and it does
# not program an FPGA.  hdl/ecc is a bare datapath: it has no AWS shell
# wrapper, no AXI interfaces, no DMA, so there is nothing for the F2 shell
# to talk to.  Producing a loadable AFI needs the aws-fpga HDK custom-logic
# wrapper first; see README.md in this directory for what that involves.
# Synthesis is the part that answers the question the cost model asks, and
# it needs no shell at all.
#
# THIS COSTS MONEY.  It launches an on-demand EC2 instance and terminates it
# when finished, including on error or interrupt.  A z1d.2xlarge is around
# $0.75/hour on demand; a run of this design should take well under an hour.
# Nothing happens until you confirm, and --dry-run shows the plan without
# spending anything.
#
# Prerequisites:
#   - awscli v2, configured (aws sts get-caller-identity must work)
#   - an EC2 key pair whose private key you have locally
#   - a subnet with a public IP and a security group allowing SSH from you
#   - a one-time subscription to the AWS FPGA Developer AMI in Marketplace
#
# Usage:
#   ./run_aws_synthesis.sh --key-name mykey --key-file ~/.ssh/mykey.pem \
#                          [--subnet subnet-xxx] [--sg sg-xxx] \
#                          [--instance-type z1d.2xlarge] [--period 3.333] \
#                          [--part xcvu47p-fsvh2892-2-e] [--region us-east-1] \
#                          [--keep] [--dry-run]

set -euo pipefail

KEY_NAME=""
KEY_FILE=""
SUBNET=""
SG=""
INSTANCE_TYPE="z1d.2xlarge"
REGION="${AWS_REGION:-us-east-1}"
PART="xcvu47p-fsvh2892-2-e"
PERIOD="3.333"          # 300 MHz, the cost model's estimate; sweep it
DRY_RUN=0
KEEP=0
OUTDIR="$(cd "$(dirname "$0")" && pwd)/reports"
HDLDIR="$(cd "$(dirname "$0")/.." && pwd)"
INSTANCE_ID=""

die() { echo "error: $*" >&2; exit 1; }
log() { printf '\033[1m==\033[0m %s\n' "$*"; }

while [ $# -gt 0 ]; do
    case "$1" in
        --key-name)      KEY_NAME="$2"; shift 2 ;;
        --key-file)      KEY_FILE="$2"; shift 2 ;;
        --subnet)        SUBNET="$2"; shift 2 ;;
        --sg)            SG="$2"; shift 2 ;;
        --instance-type) INSTANCE_TYPE="$2"; shift 2 ;;
        --region)        REGION="$2"; shift 2 ;;
        --part)          PART="$2"; shift 2 ;;
        --period)        PERIOD="$2"; shift 2 ;;
        --out)           OUTDIR="$2"; shift 2 ;;
        --keep)          KEEP=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       sed -n '2,40p' "$0"; exit 0 ;;
        *)               die "unknown option $1" ;;
    esac
done

# ---------------------------------------------------------------- #
# Terminate on every exit path, including Ctrl-C and errors.  An EC2
# instance left running because a script died is the expensive failure mode
# here, so this is armed the moment an instance exists.
# ---------------------------------------------------------------- #
cleanup() {
    local rc=$?
    if [ -n "$INSTANCE_ID" ]; then
        if [ "$KEEP" = "1" ]; then
            echo
            log "--keep given: instance $INSTANCE_ID is STILL RUNNING and still billing."
            log "terminate it with: aws ec2 terminate-instances --region $REGION --instance-ids $INSTANCE_ID"
        else
            log "terminating $INSTANCE_ID"
            aws ec2 terminate-instances --region "$REGION" \
                --instance-ids "$INSTANCE_ID" --output text >/dev/null 2>&1 \
                || echo "WARNING: could not terminate $INSTANCE_ID -- terminate it by hand" >&2
        fi
    fi
    exit $rc
}
trap cleanup EXIT INT TERM

command -v aws >/dev/null || die "awscli not found"
[ -n "$KEY_NAME" ] || die "--key-name is required"
[ -n "$KEY_FILE" ] || die "--key-file is required"
[ -f "$KEY_FILE" ] || die "key file $KEY_FILE not found"

log "checking credentials"
aws sts get-caller-identity --region "$REGION" --output text >/dev/null \
    || die "aws sts get-caller-identity failed; configure the CLI first"

log "looking up the FPGA Developer AMI"
AMI=$(aws ec2 describe-images --region "$REGION" --owners aws-marketplace \
        --filters "Name=name,Values=FPGA Developer AMI*" "Name=state,Values=available" \
        --query 'reverse(sort_by(Images,&CreationDate))[0].ImageId' --output text 2>/dev/null || true)
if [ -z "$AMI" ] || [ "$AMI" = "None" ]; then
    die "no FPGA Developer AMI found in $REGION.
   Subscribe once (free) in the AWS Marketplace, then re-run.  The AMI ships
   Vivado already licensed for AWS FPGA parts, which is why we use it rather
   than installing Vivado ourselves."
fi

cat <<PLAN

  region          $REGION
  instance type   $INSTANCE_TYPE
  AMI             $AMI  (FPGA Developer AMI)
  part            $PART
  target period   $PERIOD ns
  reports to      $OUTDIR

  This launches an on-demand instance and terminates it when done.
  A z1d.2xlarge is roughly \$0.75/hour on demand in us-east-1.

PLAN

if [ "$DRY_RUN" = "1" ]; then
    log "dry run: nothing launched"
    trap - EXIT INT TERM
    exit 0
fi

printf 'Launch it? [y/N] '
read -r reply
case "$reply" in [yY]*) ;; *) log "aborted"; trap - EXIT INT TERM; exit 0 ;; esac

RUN_ARGS=(--region "$REGION" --image-id "$AMI" --instance-type "$INSTANCE_TYPE"
          --key-name "$KEY_NAME" --count 1
          --block-device-mappings 'DeviceName=/dev/sda1,Ebs={VolumeSize=120,VolumeType=gp3,DeleteOnTermination=true}'
          --tag-specifications 'ResourceType=instance,Tags=[{Key=Name,Value=crypto-hdl-ecc-synthesis}]'
          --query 'Instances[0].InstanceId' --output text)
[ -n "$SUBNET" ] && RUN_ARGS+=(--subnet-id "$SUBNET" --associate-public-ip-address)
[ -n "$SG" ] && RUN_ARGS+=(--security-group-ids "$SG")

log "launching"
INSTANCE_ID=$(aws ec2 run-instances "${RUN_ARGS[@]}")
log "instance $INSTANCE_ID"

log "waiting for it to run"
aws ec2 wait instance-running --region "$REGION" --instance-ids "$INSTANCE_ID"
HOST=$(aws ec2 describe-instances --region "$REGION" --instance-ids "$INSTANCE_ID" \
        --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
[ "$HOST" != "None" ] || die "instance has no public IP; pass --subnet for a public subnet"
log "public ip $HOST"

SSH="ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 -i $KEY_FILE ec2-user@$HOST"

log "waiting for sshd"
for i in $(seq 1 60); do
    $SSH true 2>/dev/null && break
    [ "$i" = 60 ] && die "ssh never came up"
    sleep 10
done

log "copying rtl"
$SSH 'rm -rf ~/hdl_ecc && mkdir -p ~/hdl_ecc/aws'
scp -q -o StrictHostKeyChecking=accept-new -i "$KEY_FILE" \
    "$HDLDIR"/*.vhd "ec2-user@$HOST:~/hdl_ecc/"
scp -q -i "$KEY_FILE" "$HDLDIR/aws/synth.tcl" "ec2-user@$HOST:~/hdl_ecc/aws/"

log "running synthesis (this is the slow part)"
$SSH "bash -lc '
    set -e
    source /opt/Xilinx/Vivado/*/settings64.sh 2>/dev/null || \
      source \$(ls -d /opt/Xilinx/Vivado/*/settings64.sh | tail -1)
    cd ~/hdl_ecc
    vivado -mode batch -nojournal -log vivado.log \
           -source aws/synth.tcl -tclargs $PART $PERIOD ~/hdl_ecc/reports
'" || die "synthesis failed; re-run with --keep and ssh in to read ~/hdl_ecc/vivado.log"

log "retrieving reports"
mkdir -p "$OUTDIR"
scp -q -i "$KEY_FILE" "ec2-user@$HOST:~/hdl_ecc/reports/*" "$OUTDIR/" || true
scp -q -i "$KEY_FILE" "ec2-user@$HOST:~/hdl_ecc/vivado.log" "$OUTDIR/" || true

echo
log "results"
[ -f "$OUTDIR/summary.txt" ] && cat "$OUTDIR/summary.txt"
echo
log "full reports in $OUTDIR"
log "compare against the estimates in docs/ecc_fpga_cost_model.md and update them"
