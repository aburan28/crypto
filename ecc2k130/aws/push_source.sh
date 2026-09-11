#!/usr/bin/env bash
#
# Upload the ecc2k130 source tree to the campaign bucket so a build host needs
# no git access.  Uses the working tree (uncommitted aws/ files included) and
# excludes binaries and build products; the tarball is named by the git HEAD
# plus a content hash so build.sh can record exactly what it compiled.
#
#   ./push_source.sh            -> s3://$BUCKET/src/ecc2k130-<sha>.tar.gz, prints the key

set -euo pipefail
cd "$(dirname "$0")/.."

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
BUCKET=${BUCKET:-$STACK-$ACCOUNT}

tmp=$(mktemp -d "${TMPDIR:-/tmp}/ecc-src.XXXXXX")
tar --exclude='./ecc2k130' --exclude='./ecc2k130-cpu' --exclude='./build' --exclude='__pycache__' \
    --exclude='*.pyc' --exclude='*.ptx' --exclude='*.cubin' --exclude='*.o' \
    -czf "$tmp/src.tar.gz" .
head=$(GIT_OPTIONAL_LOCKS=0 git rev-parse --short=12 HEAD 2>/dev/null || echo nogit)
sum=$(shasum -a 256 "$tmp/src.tar.gz" | cut -c1-12)
key="src/ecc2k130-$head-$sum.tar.gz"
aws s3 cp "$tmp/src.tar.gz" "s3://$BUCKET/$key" --only-show-errors
# Record it in the campaign so an instance with no published binary can
# build from exactly this tree (bootstrap.sh), without anyone logging in.
if aws s3 cp "s3://$BUCKET/campaign.json" "$tmp/campaign.json" --only-show-errors 2>/dev/null; then
    python3 - "$tmp/campaign.json" "$key" <<'EOF'
import json, sys
path, key = sys.argv[1:]
c = json.load(open(path))
c["sourceKey"] = key
json.dump(c, open(path, "w"), indent=1, sort_keys=True)
EOF
    aws s3 cp "$tmp/campaign.json" "s3://$BUCKET/campaign.json" --only-show-errors
    echo "campaign.json sourceKey = $key" >&2
else
    echo "no campaign.json in the bucket yet (run infra.sh first); sourceKey not recorded" >&2
fi
rm -rf "$tmp"
echo "$key"
