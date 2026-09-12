# build_afi_instance.sh -- the build, as it runs on the build instance.
#
# build_afi.sh prepends a shebang and the variables (BUCKET, TAG, REGION,
# NENG, ID_W, DP_WEIGHT) and hands the result to EC2 as user data.  From
# there: the F2 HDK is cloned and set up (this downloads the shell
# checkpoint), the sources from S3 are unpacked, the CL is built with the
# HDK's own flow, the DCP tarball is uploaded, create-fpga-image is called,
# and the instance shuts down (which terminates it).  The log goes to S3
# every five minutes, so ./build_afi.sh status can follow along.

set -uo pipefail
exec > >(tee -a /var/log/ecc2k130-afi.log) 2>&1
echo "afi build $TAG starting $(date -u)"
export AWS_DEFAULT_REGION=$REGION
export HOME=${HOME:-/root}

LOGKEY="fpga/builds/$TAG/build.log"
ship() { aws s3 cp /var/log/ecc2k130-afi.log "s3://$BUCKET/$LOGKEY" --only-show-errors 2>/dev/null; }
( while sleep 300; do ship; done ) &
SHIPPER=$!
finish() {
    echo "afi build $TAG finished $(date -u)"
    kill $SHIPPER 2>/dev/null
    # Reports and the vivado log are worth having even when the build
    # failed; the utilisation of a design that did not fit is the point.
    if [ -n "${CL_DIR:-}" ] && [ -d "$CL_DIR/build/reports" ]; then
        aws s3 cp "$CL_DIR/build/reports" "s3://$BUCKET/fpga/builds/$TAG/reports/" --recursive --only-show-errors 2>/dev/null
        for l in "$CL_DIR"/build/scripts/*.vivado.log; do
            [ -f "$l" ] && aws s3 cp "$l" "s3://$BUCKET/fpga/builds/$TAG/vivado.log" --only-show-errors 2>/dev/null
        done
    fi
    ship
    [ "${NO_SHUTDOWN:-0}" = 1 ] || shutdown -h now
}
trap finish EXIT

fail() { echo "FAILED: $*"; exit 1; }

# ---- tools -----------------------------------------------------------------
# The FPGA Developer AMI has Vivado under /tools/Xilinx; a login shell would
# put it on PATH, user data does not.  Two layouts: Vivado/<ver> up to
# 2024.2, <ver>/Vivado from 2025.1.  settings64.sh also sets XILINX_VIVADO,
# which hdk_setup.sh wants.
if ! command -v vivado >/dev/null 2>&1; then
    for s in /tools/Xilinx/*/Vivado/settings64.sh /tools/Xilinx/Vivado/*/settings64.sh \
             /opt/Xilinx/*/Vivado/settings64.sh /opt/Xilinx/Vivado/*/settings64.sh; do
        [ -f "$s" ] || continue
        echo "sourcing $s"
        set +u; source "$s"; set -u
        break
    done
fi
if ! command -v vivado >/dev/null 2>&1; then
    v=$(find /tools /opt /usr/local -maxdepth 6 -type f -name vivado -path '*/bin/vivado' 2>/dev/null | sort | tail -1)
    [ -n "$v" ] && export PATH="$PATH:$(dirname "$v")"
fi
command -v vivado >/dev/null 2>&1 || {
    echo "looked for vivado; /tools: $(ls /tools 2>/dev/null | tr '\n' ' ') /tools/Xilinx: $(ls /tools/Xilinx 2>/dev/null | tr '\n' ' ') /opt: $(ls /opt 2>/dev/null | tr '\n' ' ')"
    fail "vivado not found; is this the FPGA Developer AMI?"
}
command -v aws >/dev/null 2>&1 || fail "aws cli not found"
command -v python3 >/dev/null 2>&1 || fail "python3 not found"
echo "vivado: $(vivado -version 2>/dev/null | head -1)"

# ---- HDK -------------------------------------------------------------------
WORK=/opt/ecc2k130-build
mkdir -p "$WORK"
cd "$WORK"
if [ ! -d aws-fpga ]; then
    git clone --depth 1 -b f2 https://github.com/aws/aws-fpga.git aws-fpga || fail "clone aws-fpga"
fi
cd aws-fpga
# hdk_setup.sh is written to be sourced in an interactive shell; it must not
# take our -u/-o pipefail with it.
set +u +o pipefail
source hdk_setup.sh > "$WORK/hdk_setup.log" 2>&1
rc=$?
set -u -o pipefail
grep -q "AWS HDK setup PASSED" "$WORK/hdk_setup.log" || { tail -30 "$WORK/hdk_setup.log"; fail "hdk_setup.sh (rc $rc)"; }
echo "HDK ready: shell $HDK_SHELL_DIR"

# ---- sources ---------------------------------------------------------------
cd "$WORK"
aws s3 cp "s3://$BUCKET/fpga/source.tar.gz" source.tar.gz --only-show-errors || fail "fetch fpga/source.tar.gz"
rm -rf src && mkdir src && tar xzf source.tar.gz -C src || fail "unpack source"
export CL_DIR="$WORK/src/hdl/ecc2k130/aws/cl_ecc2k130"
export ECC_RTL_DIR="$WORK/src/hdl/ecc2k130"
cd "$CL_DIR/build/scripts"
for f in aws_build_dcp_from_cl.py build_all.tcl build_level_1_cl.tcl; do
    ln -sf "$HDK_SHELL_DIR/build/scripts/$f" "$f"
done
mkdir -p "$CL_DIR/build/checkpoints" "$CL_DIR/build/reports"

# ---- build -----------------------------------------------------------------
export ECC_NENG=$NENG ECC_ID_W=$ID_W ECC_DP_WEIGHT=$DP_WEIGHT
echo "building cl_ecc2k130: $NENG engines x $((1 << ID_W)) walks, dp weight $DP_WEIGHT, tag $TAG"
python3 aws_build_dcp_from_cl.py --cl cl_ecc2k130 --tag "$TAG" || fail "aws_build_dcp_from_cl.py"

TARBALL="$CL_DIR/build/checkpoints/$TAG.Developer_CL.tar"
[ -f "$TARBALL" ] || { tail -60 "$TAG.vivado.log" 2>/dev/null; fail "no DCP tarball; see the vivado log"; }
if ls "$CL_DIR/build/checkpoints/"*VIOLATED* >/dev/null 2>&1; then
    echo "WARNING: timing was not met; the image will be flagged timingViolated in afi.json"
    TIMING=violated
else
    TIMING=met
fi

# ---- upload ----------------------------------------------------------------
PREFIX="fpga/builds/$TAG"
aws s3 cp "$TARBALL" "s3://$BUCKET/$PREFIX/$TAG.Developer_CL.tar" --only-show-errors || fail "upload tarball"
aws s3 cp "$CL_DIR/build/reports" "s3://$BUCKET/$PREFIX/reports/" --recursive --only-show-errors || true
aws s3 cp "$TAG.vivado.log" "s3://$BUCKET/$PREFIX/vivado.log" --only-show-errors || true
# The numbers ../../README.md's capacity estimate waits for, front and centre.
for r in "$CL_DIR"/build/reports/*utilization*.rpt; do
    [ -f "$r" ] && { echo "== $(basename "$r")"; grep -m1 -A12 "CLB LUTs\|Slice LUTs" "$r" || true; }
done
for r in "$CL_DIR"/build/reports/*timing*.rpt; do
    [ -f "$r" ] && { echo "== $(basename "$r")"; grep -m1 -A4 "WNS(ns)" "$r" || true; }
done

# ---- AFI -------------------------------------------------------------------
out=$(aws ec2 create-fpga-image --name "ecc2k130-$TAG" \
      --description "ECC2K-130 rho engine, $NENG engines x $((1 << ID_W)) walks, dp weight $DP_WEIGHT, timing $TIMING" \
      --input-storage-location "Bucket=$BUCKET,Key=$PREFIX/$TAG.Developer_CL.tar" \
      --logs-storage-location "Bucket=$BUCKET,Key=$PREFIX/afi-logs" \
      --tag-specifications "ResourceType=fpga-image,Tags=[{Key=Project,Value=ecc2k130},{Key=BuildTag,Value=$TAG}]" \
      --output json) || fail "create-fpga-image"
echo "$out"
python3 - "$out" "$TAG" "$NENG" "$ID_W" "$DP_WEIGHT" "$TIMING" > afi.json <<'EOF'
import json, sys
out, tag, neng, idw, dpw, timing = sys.argv[1:]
d = json.loads(out)
json.dump({"afi": d["FpgaImageId"], "agfi": d["FpgaImageGlobalId"], "tag": tag,
           "neng": int(neng), "idW": int(idw), "dpWeight": int(dpw),
           "walks": int(neng) << int(idw), "timing": timing}, sys.stdout, indent=1)
EOF
aws s3 cp afi.json "s3://$BUCKET/$PREFIX/afi.json" --only-show-errors || fail "upload afi.json"
cat afi.json
echo "submitted; ./build_afi.sh wait $TAG then ./build_afi.sh promote $TAG"
