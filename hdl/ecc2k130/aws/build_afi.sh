#!/usr/bin/env bash
#
# Build the ECC2K-130 FPGA image (an AFI) on an AWS build instance, from
# this machine with the AWS CLI.  Nothing here needs Vivado locally.
#
#   ./build_afi.sh push                 tar the engine + host sources to S3
#   ./build_afi.sh launch [TAG]         start a build instance that synthesises,
#                                       places, routes and submits the AFI
#                                       (hours; the instance shuts down after)
#   ./build_afi.sh status [TAG]         instance state, build log tail, AFI state
#   ./build_afi.sh wait   [TAG]         block until the AFI is available
#   ./build_afi.sh promote [TAG]        make it the image bootstrap_f2.sh loads
#   ./build_afi.sh list                 every build in the bucket
#
# Geometry of the image, all optional:
#   NENG       walker engines (default 48; about 8.4k LUTs, 12 RAMB36 +
#              2 RAMB18 and 4 URAM288 each; the VU47P has 1.30M LUTs,
#              2 016 RAMB36 and 960 URAM288: 112 is 73% of the LUTs, 128 is
#              83%; 48, 64 and 80 have run on a device)
#   ID_W       walks per engine = 2**ID_W (default 9: 512 walks fill the
#              FIFO's block RAM exactly)
#   DP_WEIGHT  distinguished-point cutoff baked into the image (default 34,
#              the challenge's; must equal campaign.json dpWeight)
#   CLK_MHZ    engine clock: 250, 300, 333 (default), 350, 375 or 400; or set
#              MMCM_MULT and MMCM_DIV directly (engine clock = 250 * MULT /
#              DIV, VCO = 250 * MULT within 800..1600).  Below 250 there is
#              no point; above what the routed design closes, the image is
#              flagged timing violated and its reports fail verification.
# Build instance:
#   BUILD_TYPE  default r6i.4xlarge (128 GB; Vivado on a VU47P wants > 64)
#   AMI         override the FPGA Developer AMI lookup (needs a Marketplace
#               subscription: search "AWS FPGA Developer AMI")
#   KEY_NAME    ssh key for the build instance (none = SSM only)
#   KEEP=1      leave the instance running after the build (default: it
#               shuts down and terminates)
# Shared with ecc2k130/aws/infra.sh, which must have run first:
#   AWS_DEFAULT_REGION (us-west-2), STACK (ecc2k130), BUCKET ($STACK-<account>)
#
# Layout in the bucket:
#   fpga/source.tar.gz                    what push uploads
#   fpga/builds/TAG/build.log             the instance's log, refreshed every 5 min
#   fpga/builds/TAG/TAG.Developer_CL.tar  the DCP tarball create-fpga-image ingests
#   fpga/builds/TAG/reports/              utilisation and timing
#   fpga/builds/TAG/afi.json              ids and geometry once submitted
#   fpga/afi.json                         the promoted image

set -euo pipefail
cd "$(dirname "$0")"

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-$STACK-$ACCOUNT}
ROLE=$STACK-fpga-build
PROFILE=$STACK-fpga-build
SG=$STACK-worker
BUILD_TYPE=${BUILD_TYPE:-r6i.4xlarge}
KEY_NAME=${KEY_NAME:-}
NENG=${NENG:-48}
ID_W=${ID_W:-9}
DP_WEIGHT=${DP_WEIGHT:-34}
CLK_MHZ=${CLK_MHZ:-333}
case $CLK_MHZ in
    250) M=4 D=4 ;;
    300) M=6 D=5 ;;
    333) M=4 D=3 ;;
    350) M=4.375 D=3.125 ;;
    375) M=6 D=4 ;;
    400) M=4 D=2.5 ;;
    *)   [ -n "${MMCM_MULT:-}" ] && [ -n "${MMCM_DIV:-}" ] \
             || { echo "CLK_MHZ=$CLK_MHZ has no table entry; set MMCM_MULT and MMCM_DIV" >&2; exit 1; } ;;
esac
MMCM_MULT=${MMCM_MULT:-$M}
MMCM_DIV=${MMCM_DIV:-$D}
CLK_MHZ=$(python3 -c "print(round(250 * $MMCM_MULT / $MMCM_DIV))")
REPO=$(cd ../../.. && pwd)

cmd=${1:-status}
TAG=${2:-}

lastTag() {
    aws s3 ls "s3://$BUCKET/fpga/builds/" | awk '{print $2}' | tr -d / | sort | tail -1
}
[ -n "$TAG" ] || TAG=$(lastTag || true)

case "$cmd" in
push)
    # Everything the CL build and the host program need, nothing else.
    tmp=$(mktemp -d)
    tar -C "$REPO" -czf "$tmp/source.tar.gz" \
        --exclude 'hdl/ecc2k130/host/ecc2k130-fpga' \
        --exclude 'hdl/ecc2k130/*.o' --exclude 'hdl/ecc2k130/*.cf' \
        --exclude 'hdl/ecc2k130/aws/cl_ecc2k130/build/checkpoints' \
        --exclude 'hdl/ecc2k130/aws/cl_ecc2k130/build/reports' \
        --exclude 'hdl/ecc2k130/aws/cl_ecc2k130/build/src_post_encryption' \
        hdl/ecc2k130 ecc2k130/include ecc2k130/generated ecc2k130/aws/worker.py
    aws s3 cp "$tmp/source.tar.gz" "s3://$BUCKET/fpga/source.tar.gz" --only-show-errors
    rm -rf "$tmp"
    echo "pushed source ($(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo untracked)) to s3://$BUCKET/fpga/source.tar.gz"
    ;;

launch)
    TAG=${2:-$(date -u +%Y%m%d-%H%M%S)-n$NENG-c$CLK_MHZ}
    aws s3api head-object --bucket "$BUCKET" --key fpga/source.tar.gz >/dev/null 2>&1 \
        || { echo "no source in the bucket; run ./build_afi.sh push first" >&2; exit 1; }

    # ---- role: S3 on the bucket, and the two FPGA image calls ----------
    # Without IAM rights (or with NO_ROLE=1) the instance gets a 36 h session
    # token of the calling user in its user data instead.  That is the
    # caller's own permissions on the caller's own instance for the length
    # of one build, and it is what a role would have given it, scoped wider;
    # the role is preferred when it already exists or can be made.
    ensureRole() {
        # An existing profile is enough: SSO/assumed-role callers can attach
        # it but usually cannot put-role-policy / attach-role-policy.
        aws iam get-instance-profile --instance-profile-name "$PROFILE" >/dev/null 2>&1 && return 0
        aws iam get-role --role-name "$ROLE" >/dev/null 2>&1 \
        || aws iam create-role --role-name "$ROLE" --assume-role-policy-document '{
             "Version": "2012-10-17",
             "Statement": [{"Effect": "Allow", "Principal": {"Service": "ec2.amazonaws.com"}, "Action": "sts:AssumeRole"}]
           }' >/dev/null || return 1
        aws iam put-role-policy --role-name "$ROLE" --policy-name fpga-build --policy-document "{
          \"Version\": \"2012-10-17\",
          \"Statement\": [
            {\"Effect\": \"Allow\", \"Action\": [\"s3:ListBucket\"], \"Resource\": \"arn:aws:s3:::$BUCKET\"},
            {\"Effect\": \"Allow\", \"Action\": [\"s3:GetObject\", \"s3:PutObject\"], \"Resource\": \"arn:aws:s3:::$BUCKET/*\"},
            {\"Effect\": \"Allow\", \"Action\": [\"ec2:CreateFpgaImage\", \"ec2:DescribeFpgaImages\", \"ec2:CreateTags\"], \"Resource\": \"*\"}
          ]
        }" || return 1
        aws iam attach-role-policy --role-name "$ROLE" --policy-arn arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore || return 1
        aws iam create-instance-profile --instance-profile-name "$PROFILE" >/dev/null || return 1
        aws iam add-role-to-instance-profile --instance-profile-name "$PROFILE" --role-name "$ROLE" || return 1
        echo "created instance profile $PROFILE; waiting for IAM to propagate"
        sleep 15
    }
    profileOpt=()
    credLine=""
    if [ "${NO_ROLE:-0}" != 1 ] && ensureRole; then
        profileOpt=(--iam-instance-profile "Name=$PROFILE")
        echo "instance profile $PROFILE"
    else
        echo "no IAM rights to make role $ROLE; the instance will use a 36 h session token of $(aws sts get-caller-identity --query Arn --output text)" >&2
        read -r AK SK ST <<<"$(aws sts get-session-token --duration-seconds 129600 \
                              --query 'Credentials.[AccessKeyId,SecretAccessKey,SessionToken]' --output text)"
        [ -n "$ST" ] || { echo "sts get-session-token failed" >&2; exit 1; }
        credLine="export AWS_ACCESS_KEY_ID=$AK AWS_SECRET_ACCESS_KEY=$SK AWS_SESSION_TOKEN=$ST"
    fi

    # ---- AMI: the FPGA Developer AMI carries Vivado and its licence ----
    if [ -z "${AMI:-}" ]; then
        AMI=$(aws ec2 describe-images --owners aws-marketplace \
              --filters "Name=name,Values=FPGA Developer AMI (Ubuntu)*" "Name=state,Values=available" \
                        "Name=architecture,Values=x86_64" \
              --query 'sort_by(Images,&CreationDate)[-1].ImageId' --output text)
        if [ -z "$AMI" ] || [ "$AMI" = None ]; then
            echo "FPGA Developer AMI not visible: subscribe in the AWS Marketplace or pass AMI=ami-..." >&2
            exit 1
        fi
    fi
    read -r ROOTDEV AMINAME <<<"$(aws ec2 describe-images --image-ids "$AMI" --query 'Images[0].[RootDeviceName,Name]' --output text)"
    echo "AMI $AMI ($AMINAME, root $ROOTDEV)"

    VPC=$(aws ec2 describe-vpcs --filters Name=is-default,Values=true --query 'Vpcs[0].VpcId' --output text)
    SGID=$(aws ec2 describe-security-groups --filters "Name=group-name,Values=$SG" "Name=vpc-id,Values=$VPC" \
           --query 'SecurityGroups[0].GroupId' --output text)
    [ -n "$SGID" ] && [ "$SGID" != None ] || { echo "security group $SG missing; run ecc2k130/aws/infra.sh first" >&2; exit 1; }

    # ---- user data: the whole build, unattended ------------------------
    ud=$(mktemp)
    {
        echo '#!/bin/bash'
        echo "BUCKET=$BUCKET; TAG=$TAG; REGION=$AWS_DEFAULT_REGION"
        echo "NENG=$NENG; ID_W=$ID_W; DP_WEIGHT=$DP_WEIGHT; NO_SHUTDOWN=${KEEP:-0}"
        echo "MMCM_MULT=$MMCM_MULT; MMCM_DIV=$MMCM_DIV; CLK_MHZ=$CLK_MHZ"
        [ -n "$credLine" ] && echo "$credLine"
        cat build_afi_instance.sh
    } > "$ud"
    keyOpt=()
    [ -n "$KEY_NAME" ] && keyOpt=(--key-name "$KEY_NAME")
    IID=$(aws ec2 run-instances --image-id "$AMI" --instance-type "$BUILD_TYPE" \
          ${profileOpt[@]+"${profileOpt[@]}"} --security-group-ids "$SGID" \
          --user-data "file://$ud" ${keyOpt[@]+"${keyOpt[@]}"} \
          --block-device-mappings "[{\"DeviceName\":\"$ROOTDEV\",\"Ebs\":{\"VolumeSize\":200,\"VolumeType\":\"gp3\",\"DeleteOnTermination\":true}}]" \
          --instance-initiated-shutdown-behavior terminate \
          --metadata-options HttpTokens=required,HttpPutResponseHopLimit=2 \
          --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$STACK-fpga-build-$TAG},{Key=Project,Value=$STACK},{Key=BuildTag,Value=$TAG}]" \
          --query 'Instances[0].InstanceId' --output text)
    rm -f "$ud"
    echo "build $TAG: instance $IID ($BUILD_TYPE), $NENG engines x $((1 << ID_W)) walks, dp weight $DP_WEIGHT, engine clock $CLK_MHZ MHz (MMCM $MMCM_MULT / $MMCM_DIV)"
    echo "follow with: ./build_afi.sh status $TAG   (the instance terminates itself when done)"
    ;;

status)
    [ -n "$TAG" ] || { echo "no builds yet"; exit 0; }
    echo "== build $TAG"
    aws ec2 describe-instances --filters "Name=tag:BuildTag,Values=$TAG" \
        --query 'Reservations[].Instances[].[InstanceId,State.Name,InstanceType,LaunchTime]' --output text \
        | sed 's/^/instance  /'
    if aws s3 cp "s3://$BUCKET/fpga/builds/$TAG/build.log" - 2>/dev/null | tail -n 12 | sed 's/^/log       /'; then :; fi
    if aws s3 cp "s3://$BUCKET/fpga/builds/$TAG/afi.json" - >"${TMPDIR:-/tmp}/ecc-afi.json" 2>/dev/null; then
        afi=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["afi"])' "${TMPDIR:-/tmp}/ecc-afi.json")
        aws ec2 describe-fpga-images --fpga-image-ids "$afi" \
            --query 'FpgaImages[0].[FpgaImageId,FpgaImageGlobalId,State.Code,State.Message]' --output text \
            | sed 's/^/afi       /'
    else
        echo "afi       not submitted yet"
    fi
    ;;

wait)
    [ -n "$TAG" ] || { echo "no builds yet" >&2; exit 1; }
    echo "waiting for build $TAG (Vivado takes hours; AFI generation about an hour after that)"
    while ! aws s3 cp "s3://$BUCKET/fpga/builds/$TAG/afi.json" "${TMPDIR:-/tmp}/ecc-afi.json" --only-show-errors 2>/dev/null; do
        if [ "$(aws ec2 describe-instances --filters "Name=tag:BuildTag,Values=$TAG" "Name=instance-state-name,Values=pending,running,stopping" \
                 --query 'length(Reservations[].Instances[])' --output text)" = 0 ]; then
            echo "build instance is gone and no AFI was submitted; read fpga/builds/$TAG/build.log" >&2
            exit 1
        fi
        sleep 300
    done
    afi=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["afi"])' "${TMPDIR:-/tmp}/ecc-afi.json")
    while :; do
        state=$(aws ec2 describe-fpga-images --fpga-image-ids "$afi" --query 'FpgaImages[0].State.Code' --output text)
        case "$state" in
            available) echo "AFI $afi available"; exit 0 ;;
            failed) aws ec2 describe-fpga-images --fpga-image-ids "$afi" --query 'FpgaImages[0].State.Message' --output text >&2; exit 1 ;;
            *) echo "  $(date -u +%H:%M) $state"; sleep 120 ;;
        esac
    done
    ;;

promote)
    [ -n "$TAG" ] || { echo "no builds yet" >&2; exit 1; }
    aws s3 cp "s3://$BUCKET/fpga/builds/$TAG/afi.json" "${TMPDIR:-/tmp}/ecc-afi.json" --only-show-errors
    afi=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["afi"])' "${TMPDIR:-/tmp}/ecc-afi.json")
    state=$(aws ec2 describe-fpga-images --fpga-image-ids "$afi" --query 'FpgaImages[0].State.Code' --output text)
    [ "$state" = available ] || { echo "AFI $afi is $state, not available" >&2; exit 1; }
    aws s3 cp "${TMPDIR:-/tmp}/ecc-afi.json" "s3://$BUCKET/fpga/afi.json" --only-show-errors
    python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print("promoted %s: %s, %d engines x %d walks, dp weight %d" % (d["tag"], d["agfi"], d["neng"], 1 << d["idW"], d["dpWeight"]))' "${TMPDIR:-/tmp}/ecc-afi.json"
    echo "F2 instances started from now on load it; running ones reload on their next restart"
    ;;

list)
    for t in $(aws s3 ls "s3://$BUCKET/fpga/builds/" | awk '{print $2}' | tr -d /); do
        if aws s3 cp "s3://$BUCKET/fpga/builds/$t/afi.json" - 2>/dev/null \
            | python3 -c 'import json,sys; d=json.load(sys.stdin); print("%-28s %s  %3d eng x %4d walks  w<=%d" % (d["tag"], d["agfi"], d["neng"], 1 << d["idW"], d["dpWeight"]))'; then :
        else
            echo "$t  (not submitted)"
        fi
    done
    if aws s3 cp "s3://$BUCKET/fpga/afi.json" - 2>/dev/null | python3 -c 'import json,sys; print("promoted: " + json.load(sys.stdin)["tag"])'; then :; fi
    ;;

*)
    sed -n '3,15p' "$0"; exit 1 ;;
esac
