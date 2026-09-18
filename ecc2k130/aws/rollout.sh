#!/usr/bin/env bash
#
# Stage / activate / wait for a CUDA-kernel rollout. A compile is not a
# fleet flip: build.sh --stage publishes bin/<sha>/ and this script is
# what rewrites campaign.json. Activate refuses frozen geometry (see
# rollout.py). Workers poll binaryKey on the heartbeat and SIGTERM /
# checkpoint / reload themselves. No SSM, no TerminateInstances.
#
#   ./rollout.sh stage <prefix>              mark an already-published prefix staged
#   ./rollout.sh status
#   ./rollout.sh versions                    list immutable kernels/<n>.json
#   ./rollout.sh activate <prefix>           next kernelVersion, then point
#   ./rollout.sh activate --from-version N   re-activate that record as a new version
#   ./rollout.sh wait [prefix]               until walking slots heartbeat the new key
#
# Prefix is bin/<16-hex> (no trailing /ecc2k130).
# kernelProtocol ecc2k-kernel-v1 is the client pointer. It does not
# migrate the DP corpus and is not storageProtocol.
set -euo pipefail
cd "$(dirname "$0")"

export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
STACK=${STACK:-ecc2k130}
if [ -z "${BUCKET:-}" ]; then
    ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
    BUCKET=$STACK-$ACCOUNT
fi

cmd=${1:-}
shift || true
HERE=$(pwd)

need_prefix() {
    local p=${1:-}
    [ -n "$p" ] || { echo "usage: $0 $cmd bin/<sha>" >&2; exit 2; }
    p=${p#s3://$BUCKET/}
    p=${p%/ecc2k130}
    p=${p%/}
    echo "$p"
}

s3get() {  # key dest
    aws s3 cp "s3://$BUCKET/$1" "$2" --only-show-errors
}

prefix_complete() {
    local prefix=$1 dest=$2
    local f
    for f in ecc2k130 ecc2k130-cpu test-packed-cuda test-packed-storage-cuda \
             test-shared-sigma-cuda manifest.json libgomp.so.1; do
        aws s3api head-object --bucket "$BUCKET" --key "$prefix/$f" >/dev/null \
            || { echo "prefix $prefix is missing $f" >&2; return 1; }
    done
    s3get "$prefix/manifest.json" "$dest"
}

write_record() {
    local prefix=$1 dest=$2 status=${3:-staged}
    python3 - "$HERE/rollout.py" "$dest" "$prefix" "$status" <<'EOF'
import json, sys, importlib.util
spec = importlib.util.spec_from_file_location("rollout", sys.argv[1])
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
dest, prefix, status = sys.argv[2:]
man = json.load(open(dest))
rec = mod.stagedRecord(prefix, man, status)
json.dump(rec, open(dest, "w"), indent=1)
print(rec["binaryKey"])
EOF
}

cmd_stage() {
    local prefix dest
    prefix=$(need_prefix "${1:-}")
    dest=$(mktemp)
    prefix_complete "$prefix" "$dest" || { rm -f "$dest"; exit 1; }
    write_record "$prefix" "$dest" staged >/dev/null
    aws s3 cp "$dest" "s3://$BUCKET/rollouts/$(basename "$prefix").json" --only-show-errors
    aws s3 cp "$dest" "s3://$BUCKET/rollouts/current.json" --only-show-errors
    echo "staged s3://$BUCKET/$prefix/ as rollouts/$(basename "$prefix").json"
    echo "activate with: $0 activate $prefix"
    rm -f "$dest"
}

cmd_status() {
    local work live live_key staged
    work=$(mktemp -d)
    s3get campaign.json "$work/campaign.json"
    python3 - "$work/campaign.json" <<'EOF'
import json, sys
c = json.load(open(sys.argv[1]))
print("live  %s  %s kernelVersion=%s" % (
    c.get("binaryKey", ""), c.get("kernelProtocol") or "(unversioned)",
    c.get("kernelVersion") if c.get("kernelVersion") not in (None, "") else 0))
EOF
    live_key=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("binaryKey",""))' "$work/campaign.json")
    if aws s3api head-object --bucket "$BUCKET" --key rollouts/current.json >/dev/null 2>&1; then
        s3get rollouts/current.json "$work/current.json"
        python3 - "$work/current.json" <<'EOF'
import json, sys
r = json.load(open(sys.argv[1]))
print("stage %s  %s  arches=%s  status=%s" % (
    r.get("prefix",""), r.get("binaryKey",""), " ".join(r.get("arches") or []), r.get("status","")))
EOF
    else
        echo "stage (none)"
    fi
    if [ -n "$live_key" ]; then
        aws s3api list-objects-v2 --bucket "$BUCKET" --prefix slots/ \
            --query 'Contents[].Key' --output json > "$work/keys.json" || true
        python3 - "$work/keys.json" "$BUCKET" "$live_key" "$HERE/rollout.py" <<'EOF'
import json, os, subprocess, sys, tempfile, importlib.util
keys, bucket, live, rollout = sys.argv[1:5]
spec = importlib.util.spec_from_file_location("rollout", rollout)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
names = json.load(open(keys)) or []
slots = []
for key in names:
    name = os.path.basename(key)
    if not (name.startswith("slot-") and name.endswith(".json")):
        continue
    dest = tempfile.mktemp()
    r = subprocess.run(["aws", "s3", "cp", "s3://%s/%s" % (bucket, key), dest, "--only-show-errors"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        continue
    it = json.load(open(dest))
    it["slot"] = int(name[5:-5])
    slots.append(it)
    os.remove(dest)
rep = mod.walkingAdoption(slots, live)
print("walking %d  adopted %d  stale %d  unknown %d" % (
    rep["walking"], rep["adopted"], rep["stale"], rep["unknown"]))
EOF
    fi
    rm -rf "$work"
}

cmd_versions() {
    local work key name
    work=$(mktemp -d)
    aws s3api list-objects-v2 --bucket "$BUCKET" --prefix kernels/ \
        --query 'Contents[].Key' --output json > "$work/keys.json" || echo '[]' > "$work/keys.json"
    python3 - "$work/keys.json" "$BUCKET" <<'EOF'
import json, os, subprocess, sys
names = json.load(open(sys.argv[1])) or []
bucket = sys.argv[2]
rows = []
for key in names:
    name = os.path.basename(key)
    if not name.endswith(".json") or name == "current.json":
        continue
    dest = os.path.join("/tmp", "ecc-k-" + name)
    r = subprocess.run(["aws", "s3", "cp", "s3://%s/%s" % (bucket, key), dest, "--only-show-errors"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        continue
    rec = json.load(open(dest))
    rows.append(rec)
    os.remove(dest)
rows.sort(key=lambda r: int(r.get("kernelVersion") or 0))
if not rows:
    print("kernels (none)")
for rec in rows:
    print("v%-4s %s  %s  %s" % (
        rec.get("kernelVersion"), rec.get("status", ""), rec.get("prefix", ""),
        rec.get("binarySha256", "")[:16]))
EOF
    rm -rf "$work"
}

cmd_activate() {
    local prefix work live_key live_prefix from_ver extra allow
    from_ver=""
    if [ "${1:-}" = "--from-version" ]; then
        from_ver=${2:?usage: $0 activate --from-version N}
        shift 2 || true
        work=$(mktemp -d)
        s3get "kernels/$from_ver.json" "$work/from.json"
        prefix=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["prefix"])' "$work/from.json")
        rm -rf "$work"
        echo "re-activating kernelVersion $from_ver prefix $prefix as a new version"
    fi
    prefix=$(need_prefix "${prefix:-${1:-}}")
    work=$(mktemp -d)
    prefix_complete "$prefix" "$work/staged.json" || { rm -rf "$work"; exit 1; }
    s3get campaign.json "$work/campaign.json"
    live_key=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("binaryKey",""))' "$work/campaign.json")
    live_prefix=${live_key%/ecc2k130}
    extra=()
    if [ -n "${KERNEL_ALLOW_STRICT:-}" ]; then
        extra+=(--allow-strict)
    fi
    if [ -n "$live_prefix" ] && aws s3api head-object --bucket "$BUCKET" --key "$live_prefix/manifest.json" >/dev/null 2>&1; then
        s3get "$live_prefix/manifest.json" "$work/live-manifest.json"
        extra+=(--live-manifest "$work/live-manifest.json")
    fi
    python3 "$HERE/rollout.py" check --campaign "$work/campaign.json" --staged "$work/staged.json" "${extra[@]}" \
        || { echo "activate refused; campaign.json not rewritten" >&2; rm -rf "$work"; exit 2; }
    kver=$(python3 "$HERE/rollout.py" next-version --campaign "$work/campaign.json" | python3 -c 'import json,sys; print(json.load(sys.stdin)["nextKernelVersion"])')
    python3 "$HERE/rollout.py" apply --campaign "$work/campaign.json" --staged "$work/staged.json" \
        --prefix "$prefix" --kernel-version "$kver" --record "$work/kernel.json" "${extra[@]}"
    if ! aws s3api put-object --bucket "$BUCKET" --key "kernels/$kver.json" \
            --body "$work/kernel.json" --content-type application/json \
            --if-none-match "*" >/dev/null; then
        echo "activate refused: kernels/$kver.json already exists; campaign.json not rewritten" >&2
        rm -rf "$work"
        exit 1
    fi
    aws s3 cp "$work/kernel.json" "s3://$BUCKET/kernels/current.json" --only-show-errors
    aws s3 cp "$work/kernel.json" "s3://$BUCKET/rollouts/$(basename "$prefix").json" --only-show-errors
    aws s3 cp "$work/kernel.json" "s3://$BUCKET/rollouts/current.json" --only-show-errors
    aws s3 cp "$work/campaign.json" "s3://$BUCKET/campaign.json" --only-show-errors
    echo "campaign.json now selects $prefix/ecc2k130 (ecc2k-kernel-v1 kernelVersion=$kver)"
    echo "workers will SIGTERM, checkpoint and reload on the next heartbeat"
    echo "watch with: $0 wait $prefix"
    rm -rf "$work"
}

cmd_wait() {
    local prefix live_key want_key want_sha work deadline now
    if [ -n "${1:-}" ]; then
        prefix=$(need_prefix "$1")
        want_key="$prefix/ecc2k130"
    else
        work=$(mktemp -d)
        s3get campaign.json "$work/campaign.json"
        want_key=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("binaryKey",""))' "$work/campaign.json")
        rm -rf "$work"
        prefix=${want_key%/ecc2k130}
    fi
    [ -n "$want_key" ] || { echo "no binaryKey to wait for" >&2; exit 2; }
    want_sha=""
    want_ver=""
    if aws s3api head-object --bucket "$BUCKET" --key campaign.json >/dev/null 2>&1; then
        work=$(mktemp)
        s3get campaign.json "$work"
        want_ver=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("kernelVersion") or "")' "$work")
        rm -f "$work"
    fi
    if aws s3api head-object --bucket "$BUCKET" --key "$prefix/manifest.json" >/dev/null 2>&1; then
        work=$(mktemp)
        s3get "$prefix/manifest.json" "$work"
        want_sha=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("binarySha256",""))' "$work")
        rm -f "$work"
    fi
    deadline=$(( $(date +%s) + ${WAIT_SECONDS:-3600} ))
    echo "waiting for walking slots to heartbeat $want_key"
    while :; do
        work=$(mktemp -d)
        aws s3api list-objects-v2 --bucket "$BUCKET" --prefix slots/ \
            --query 'Contents[].Key' --output json > "$work/keys.json" || echo '[]' > "$work/keys.json"
        python3 - "$work/keys.json" "$BUCKET" "$want_key" "$want_sha" "$want_ver" "$HERE/rollout.py" > "$work/rep.json" <<'EOF'
import json, os, subprocess, sys, tempfile, importlib.util
keys, bucket, want, sha, ver, rollout = sys.argv[1:7]
spec = importlib.util.spec_from_file_location("rollout", rollout)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
names = json.load(open(keys)) or []
slots = []
for key in names:
    name = os.path.basename(key)
    if not (name.startswith("slot-") and name.endswith(".json")):
        continue
    dest = tempfile.mktemp()
    r = subprocess.run(["aws", "s3", "cp", "s3://%s/%s" % (bucket, key), dest, "--only-show-errors"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        continue
    it = json.load(open(dest))
    it["slot"] = int(name[5:-5])
    slots.append(it)
    os.remove(dest)
json.dump(mod.walkingAdoption(slots, want, sha, ver or None), sys.stdout)
sys.stdout.write("\n")
EOF
        python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); print("walking %(walking)s  adopted %(adopted)s  stale %(stale)s  unknown %(unknown)s"%r); sys.exit(0 if r["done"] or r["walking"]==0 else 1)' "$work/rep.json" \
            && { echo "rollout adopted (or no walking slots)"; rm -rf "$work"; return 0; }
        rm -rf "$work"
        now=$(date +%s)
        [ "$now" -lt "$deadline" ] || { echo "timed out waiting for adoption" >&2; return 1; }
        sleep "${WAIT_POLL:-30}"
    done
}

case "$cmd" in
    stage) cmd_stage "${1:-}" ;;
    status) cmd_status ;;
    versions) cmd_versions ;;
    activate) cmd_activate "$@" ;;
    wait) cmd_wait "${1:-}" ;;
    *) echo "usage: $0 {stage|status|versions|activate|wait} [bin/<sha>|--from-version N]" >&2; exit 2 ;;
esac
