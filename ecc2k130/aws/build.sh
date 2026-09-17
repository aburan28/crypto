#!/usr/bin/env bash
#
# Build the campaign binaries on a Linux host with Docker, in the same image
# the audited Modal preset used (nvidia/cuda:13.3.1-devel-ubuntu24.04), and
# publish them to the campaign bucket.  No GPU is needed to compile.
#
#   BUCKET=... ./build.sh s3-key-of-source-tarball     (from push_source.sh)
#   BUCKET=... ./build.sh /path/to/ecc2k130            (a checkout on this host)
#
# Produces s3://$BUCKET/bin/<sha>/{ecc2k130,ecc2k130-cpu,test-packed-cuda,
# test-packed-storage-cuda,test-shared-sigma-cuda,manifest.json} and points
# campaign.json at them.  The client is compiled with the RTX PRO 6000 preset
# knobs (batch 16, 256-thread blocks, minBlocks 2, native carryless products,
# 256-worker tiles, compact storage, weighted prefixes and shared Frobenius
# masks) for sm_120 only unless ARCHES says otherwise.  Native carryless
# multiplication needs CUDA 13.3 or newer, which is why the image moved off
# 13.0; the toolkit and knobs here are the ones ../RTX-PRO6000.md measured at
# 14.1 B iterations/s collecting.
#
# CLMAD is the one preset knob that does NOT travel with every ARCHES value.
# NVIDIA documents `clmad` for sm_80 and later. Blackwell pays for it on a
# measured pipe balance (87.3% ALU / 51.4% carryless on sm_120). Ada (sm_89,
# EC2 g6/g6e) now has its own receipt: CLMAD=1 is 1.811× the software product
# on an L40S and 1.881× on an L4 (../benchmarks/preblackwell/summary.json).
# Turing (sm_75, g4dn) cannot issue the instruction. Ampere (80/86) and
# Hopper (90) remain unmeasured, so they still default off. CLMAD defaults
# to 1 for Blackwell and/or Ada and to 0 the moment 75, 80, 86 or 90 is in
# ARCHES; set CLMAD=1 explicitly to override after a receipt.
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
ARCHES=${ARCHES:-120}                 # "120", "89", "75", or "75 89 120"
# Default on for Blackwell (120) and Ada (89). Off if Turing or an
# unmeasured pre-Blackwell arch is in the set; see the header comment.
# sm_75 is the EC2 g4dn (T4) client; clmad is not a T4 instruction.
case " $ARCHES " in
    *" 75 "*|*" 80 "*|*" 86 "*|*" 90 "*) clmadDefault=0 ;;
    *) clmadDefault=1 ;;
esac
CLMAD=${CLMAD:-$clmadDefault}
if [ "$CLMAD" != 0 ] && [ "$CLMAD" != 1 ]; then
    echo "CLMAD must be 0 or 1" >&2; exit 2
fi
CUDA_IMAGE=${CUDA_IMAGE:-nvidia/cuda:13.3.1-devel-ubuntu24.04}
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
# The iteration function comes from campaign.json ("walk": "sigma" or "table",
# see protocol.WALKS) so the binary and the contract that names it cannot
# disagree; the file is fetched again below to point it at the build.
aws s3 cp "s3://$BUCKET/campaign.json" "$work/campaign.json" --only-show-errors
walk=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("walk", "sigma"))' "$work/campaign.json")
case "$walk" in
    sigma) walkTable=0 ;;
    table) walkTable=1 ;;
    *) echo "campaign.json names an unknown walk: $walk" >&2; exit 1 ;;
esac
knobs="BATCH=16 THREADS=256 MINBLOCKS=2 PACKED_SINGLE_PRODUCT=1 PACKED_CACHE_DENOM=1 PACKED_BY_VALUE=1 \
PACKED_PERM_SIGMA=3 PACKED_POLY_CHAIN=1 PACKED_UNROLL_INV=1 PACKED_PAIR_PRODUCTS=1 PACKED_POLY_STATE=1 \
PACKED_DIRECT_REDUCE=1 PACKED_GENERATED_PRODUCT=1 PACKED_CLMAD=$CLMAD PACKED_STATE_TILE=256 \
PACKED_WEIGHTED_PREFIX=2 PACKED_COMPACT_STATE=1 PACKED_SHARED_SIGMA=1 WALK_TABLE=$walkTable TABLE_BRANCHES=8"
# The fixtures must be built with the arithmetic they are meant to check, so
# take their -D flags from the knobs above rather than from a second list.
defs="-DECC_STREAM_KARAT=0 -DECC_SMEM_SPILL=0"
for kv in $knobs; do defs="$defs -DECC_${kv%%=*}=${kv#*=}"; done

echo "building in $CUDA_IMAGE for sm_{$ARCHES} (ptxas takes several minutes per architecture)"
# The Deep Learning AMI's docker can reset mid-pull on a fresh box (the first
# g6 fleet died on "connection reset by peer" / "client connection is closing").
# Three attempts; the published prefix is only pointed at after a successful
# compile.
docker_ok=0
for attempt in 1 2 3; do
    echo "docker compile attempt $attempt/3"
    if docker run --rm -v "$work/src:/src" -w /src "$CUDA_IMAGE" bash -euo pipefail -c "
    apt-get update -q >/dev/null && apt-get install -y -q build-essential python3 >/dev/null
    cd codegen && python3 gen.py --out ../generated --leaf 0 && cd ..
    make gpu ARCH='$gencode' $knobs 2>&1 | tail -n 40
    make cpu MARCH=x86-64-v3 $knobs 2>&1 | tail -n 5
    mkdir -p build
    nvcc -O3 -std=c++17 $gencode $defs src/testpackedcuda.cu -o build/test-packed-cuda 2>&1 | tail -n 5
    nvcc -O3 -std=c++17 $gencode $defs src/testpackedstatecuda.cu -o build/test-packed-storage-cuda 2>&1 | tail -n 5
    nvcc -O3 -std=c++17 $gencode $defs src/testsharedsigmacuda.cu -o build/test-shared-sigma-cuda 2>&1 | tail -n 5
    if [ $walkTable = 1 ]; then
        nvcc -O3 -std=c++17 $gencode $defs -Xcompiler -fopenmp src/testtablewalkcuda.cu -o build/test-table-walk-cuda -lgomp 2>&1 | tail -n 5
    fi
    nvcc --version | tail -n 2 > build/nvcc.txt
    # The client links libgomp dynamically; ship the container's copy so the
    # instance needs no compiler packages.
    cp /usr/lib/x86_64-linux-gnu/libgomp.so.1 build/libgomp.so.1
    chown -R $(id -u):$(id -g) . 2>/dev/null || true
"
    then
        docker_ok=1
        break
    fi
    echo "docker compile attempt $attempt failed; waiting 20s"
    sleep 20
done
if [ "$docker_ok" != 1 ]; then
    echo "docker compile failed after 3 attempts" >&2
    exit 1
fi

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
# The source sha alone no longer identifies a build: ARCHES and CLMAD now vary
# independently of it, and two builds that differ only in those would publish to
# the same key and silently overwrite each other -- leaving campaign.json
# pointing at a binary whose manifest describes the other one.  Bind the key to
# everything that reaches the compiler.
buildSha=$(printf '%s\n%s\n%s\n' "$sha" "$ARCHES" "$(echo $knobs)" | sha256sum | cut -d' ' -f1)
short=${buildSha:0:16}
binSha=$(sha256sum "$work/src/ecc2k130" | cut -d' ' -f1)
hostSha=$(sha256sum "$work/src/ecc2k130-cpu" | cut -d' ' -f1)
python3 - "$work/src" "$sha" "$binSha" "$ARCHES" "$SRC" "$knobs" "$buildSha" <<'EOF' > "$work/manifest.json"
import json, sys, time
src, sha, binSha, arches, origin, knobs, buildSha = sys.argv[1:]
print(json.dumps({"sourceSha256": sha, "buildSha256": buildSha,
                  "binarySha256": binSha, "arches": arches.split(),
                  "origin": origin, "nvcc": open(src + "/build/nvcc.txt").read().strip(),
                  "knobs": " ".join(knobs.split()),
                  "builtAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1))
EOF
prefix="bin/$short"
published="ecc2k130 ecc2k130-cpu test-packed-cuda test-packed-storage-cuda test-shared-sigma-cuda libgomp.so.1"
for f in ecc2k130 ecc2k130-cpu build/test-packed-cuda build/test-packed-storage-cuda \
         build/test-shared-sigma-cuda build/libgomp.so.1; do
    [ -s "$work/src/$f" ] || { echo "build did not produce $f; not pointing campaign.json at $prefix" >&2; exit 1; }
    aws s3 cp "$work/src/$f" "s3://$BUCKET/$prefix/$(basename "$f")" --only-show-errors
done
if [ $walkTable = 1 ]; then
    [ -s "$work/src/build/test-table-walk-cuda" ] || { echo "build did not produce test-table-walk-cuda" >&2; exit 1; }
    aws s3 cp "$work/src/build/test-table-walk-cuda" "s3://$BUCKET/$prefix/test-table-walk-cuda" --only-show-errors
    published="$published test-table-walk-cuda"
fi
aws s3 cp "$work/manifest.json" "s3://$BUCKET/$prefix/manifest.json" --only-show-errors
published="$published manifest.json"
for f in $published; do
    aws s3api head-object --bucket "$BUCKET" --key "$prefix/$f" >/dev/null \
        || { echo "refusing to point campaign.json at $prefix: $f missing after upload" >&2; exit 1; }
done

# Point the campaign at this build.  Geometry in campaign.json must match the
# knobs above; they are the audited preset, so only binaryKey moves.
aws s3 cp "s3://$BUCKET/campaign.json" "$work/campaign.json" --only-show-errors
python3 - "$work/campaign.json" "$prefix" "$sha" "$knobs" "$binSha" "$hostSha" <<'EOF'
import json, sys
path, prefix, sha, knobs, binSha, hostSha = sys.argv[1:]
c = json.load(open(path))
c["binaryKey"] = prefix + "/ecc2k130"
c["hostBinaryKey"] = prefix + "/ecc2k130-cpu"
c["sourceSha256"] = sha
c["binarySha256"] = binSha
c["hostBinarySha256"] = hostSha
built = dict(kv.split("=", 1) for kv in knobs.split())
geometry = (int(built["BATCH"]), int(built["THREADS"]), int(built["MINBLOCKS"]))
assert (c["batch"], c["blockThreads"], c["minBlocks"]) == geometry, \
    "campaign.json geometry %r differs from the build knobs %r" % (
        (c["batch"], c["blockThreads"], c["minBlocks"]), geometry)
assert {"sigma": "0", "table": "1"}[c.get("walk", "sigma")] == built["WALK_TABLE"], \
    "campaign.json walk %r differs from the build's WALK_TABLE=%s" % (c.get("walk", "sigma"), built["WALK_TABLE"])
json.dump(c, open(path, "w"), indent=1, sort_keys=True)
EOF
aws s3 cp "$work/campaign.json" "s3://$BUCKET/campaign.json" --only-show-errors
echo "published s3://$BUCKET/$prefix/ (source $sha, sm_{$ARCHES}, CLMAD=$CLMAD)"
echo "campaign.json now selects $prefix/ecc2k130"
rm -rf "$work"
