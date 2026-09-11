#!/usr/bin/env bash
#
# Build the campaign binaries on a Linux host with Docker, in the same image
# the audited Modal preset used (nvidia/cuda:13.0.0-devel-ubuntu24.04), and
# publish them to the campaign bucket.  No GPU is needed to compile.
#
#   BUCKET=... ./build.sh s3-key-of-source-tarball     (from push_source.sh)
#   BUCKET=... ./build.sh /path/to/ecc2k130            (a checkout on this host)
#
# Produces s3://$BUCKET/bin/<sha>/{ecc2k130,ecc2k130-cpu,test-packed-cuda,
# manifest.json} and points campaign.json at them.  The client is compiled
# with the RTX PRO 6000 preset knobs (batch 32, 256-thread blocks, minBlocks 2,
# every packed arithmetic option) for sm_120 only unless ARCHES says otherwise.
#
# Typical use: run on the pilot g7e instance over ssh/SSM, where Docker and
# the AWS CLI are present and the instance role can write to the bucket.

set -euo pipefail
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
if [ -z "${BUCKET:-}" ]; then
    ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
    BUCKET=$STACK-$ACCOUNT
fi
ARCHES=${ARCHES:-120}                 # "120" or "89 90 120"
CUDA_IMAGE=${CUDA_IMAGE:-nvidia/cuda:13.0.0-devel-ubuntu24.04}
SRC=${1:?source: s3 key of a push_source.sh tarball, or a local ecc2k130 directory}

work=$(mktemp -d "${TMPDIR:-/tmp}/ecc-build.XXXXXX")
if [ -d "$SRC" ]; then
    cp -R "$SRC/." "$work/src"
else
    aws s3 cp "s3://$BUCKET/$SRC" "$work/src.tar.gz" --only-show-errors
    mkdir -p "$work/src"
    tar -xzf "$work/src.tar.gz" -C "$work/src"
fi
rm -rf "$work/src/build" "$work/src/ecc2k130" "$work/src/ecc2k130-cpu"

gencode=""
for a in $ARCHES; do gencode="$gencode -gencode arch=compute_$a,code=sm_$a"; done
knobs="BATCH=32 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 \
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
PACKED_DIRECT_REDUCE=1"

echo "building in $CUDA_IMAGE for sm_{$ARCHES} (ptxas takes several minutes per architecture)"
docker run --rm -v "$work/src:/src" -w /src "$CUDA_IMAGE" bash -euo pipefail -c "
    apt-get update -q >/dev/null && apt-get install -y -q build-essential python3 >/dev/null
    cd codegen && python3 gen.py --out ../generated --leaf 0 && cd ..
    make gpu ARCH='$gencode' $knobs 2>&1 | tail -n 40
    make cpu MARCH=x86-64-v3 $knobs 2>&1 | tail -n 5
    mkdir -p build
    nvcc -O3 -std=c++17 $gencode \
        -DECC_BATCH=32 -DECC_THREADS=256 -DECC_MINBLOCKS=2 -DECC_STREAM_KARAT=0 -DECC_SMEM_SPILL=0 \
        -DECC_PACKED_SINGLE_PRODUCT=1 -DECC_PACKED_CACHE_DENOM=1 -DECC_PACKED_BY_VALUE=1 \
        -DECC_PACKED_PERM_SIGMA=3 -DECC_PACKED_POLY_CHAIN=1 -DECC_PACKED_UNROLL_INV=1 \
        -DECC_PACKED_PAIR_PRODUCTS=1 -DECC_PACKED_POLY_STATE=1 -DECC_PACKED_DIRECT_REDUCE=1 \
        src/testpackedcuda.cu -o build/test-packed-cuda 2>&1 | tail -n 5
    nvcc --version | tail -n 2 > build/nvcc.txt
    # The client links libgomp dynamically; ship the container's copy so the
    # instance needs no compiler packages.
    cp /usr/lib/x86_64-linux-gnu/libgomp.so.1 build/libgomp.so.1
    chown -R $(id -u):$(id -g) . 2>/dev/null || true
"

# Source identity: the same recipe as modal_app.benchmarkIdentity, so a number
# measured on Modal and a number measured here can be tied to the same code.
sha=$(cd "$work/src" && python3 - <<'EOF'
import hashlib, pathlib
root = pathlib.Path('.')
h = hashlib.sha256()
paths = [root / 'Makefile', root / 'modal_app.py']
for folder in ('include', 'src', 'generated', 'codegen'):
    paths += sorted(p for p in (root / folder).rglob('*') if p.suffix in ('.h', '.cuh', '.cu', '.cpp', '.py'))
for p in paths:
    h.update(str(p).encode() + b'\0' + p.read_bytes() + b'\0')
print(h.hexdigest())
EOF
)
short=${sha:0:16}
binSha=$(sha256sum "$work/src/ecc2k130" | cut -d' ' -f1)
python3 - "$work/src" "$sha" "$binSha" "$ARCHES" "$SRC" <<'EOF' > "$work/manifest.json"
import json, sys, time
src, sha, binSha, arches, origin = sys.argv[1:]
print(json.dumps({"sourceSha256": sha, "binarySha256": binSha, "arches": arches.split(),
                  "origin": origin, "nvcc": open(src + "/build/nvcc.txt").read().strip(),
                  "knobs": "B32/T256/min2 single/cache/value/sigma3/poly/unroll/pair/polystate/direct",
                  "builtAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1))
EOF
prefix="bin/$short"
for f in ecc2k130 ecc2k130-cpu build/test-packed-cuda build/libgomp.so.1; do
    aws s3 cp "$work/src/$f" "s3://$BUCKET/$prefix/$(basename "$f")" --only-show-errors
done
aws s3 cp "$work/manifest.json" "s3://$BUCKET/$prefix/manifest.json" --only-show-errors

# Point the campaign at this build.  Geometry in campaign.json must match the
# knobs above; they are the audited preset, so only binaryKey moves.
aws s3 cp "s3://$BUCKET/campaign.json" "$work/campaign.json" --only-show-errors
python3 - "$work/campaign.json" "$prefix" "$sha" <<'EOF'
import json, sys
path, prefix, sha = sys.argv[1:]
c = json.load(open(path))
c["binaryKey"] = prefix + "/ecc2k130"
c["hostBinaryKey"] = prefix + "/ecc2k130-cpu"
c["sourceSha256"] = sha
assert (c["batch"], c["blockThreads"], c["minBlocks"]) == (32, 256, 2), "campaign.json geometry differs from the build knobs"
json.dump(c, open(path, "w"), indent=1, sort_keys=True)
EOF
aws s3 cp "$work/campaign.json" "s3://$BUCKET/campaign.json" --only-show-errors
echo "published s3://$BUCKET/$prefix/ (source $sha)"
echo "campaign.json now selects $prefix/ecc2k130"
rm -rf "$work"
