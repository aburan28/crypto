#!/usr/bin/env bash
#
# End-to-end test of `ephemeral` without an AWS account.
#
# tests/fake-aws stands in for the aws CLI over a scratch directory, and its
# `ec2 run-instances` boots an amazonlinux:2023 container that runs the real
# user data and runner.sh: dnf, a shallow sparse fetch of the pinned commit
# from GitHub, the local patch, the command, log streaming, output upload and
# the self-termination. Needs docker and network access; takes a few minutes.
#
#   taskq/deploy/ephemeral/tests/local_e2e.sh

set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CLI=$HERE/../ephemeral
ROOT=$(git -C "$HERE" rev-parse --show-toplevel)

FAKE_AWS=$(mktemp -d)
export FAKE_AWS XDG_CACHE_HOME=$FAKE_AWS/cache EPHEMERAL_POLL=2
mkdir -p "$FAKE_AWS/bin" "$FAKE_AWS/ssm/ephemeral"
ln -s "$HERE/fake-aws" "$FAKE_AWS/bin/aws"
export PATH=$FAKE_AWS/bin:$PATH
cleanup() {
    docker ps -aq --filter name=ephemeral-test- | xargs -r docker rm -f >/dev/null
    # Containers write as root.
    docker run --rm -v "$FAKE_AWS:/x" amazonlinux:2023 rm -rf /x/s3 /x/ec2 >/dev/null 2>&1 || true
    rm -rf "$FAKE_AWS"
}
trap cleanup EXIT

cat >"$FAKE_AWS/ssm/ephemeral/config" <<'EOF'
{"region":"us-east-1","bucket":"test-runs","launch_template":"lt-0test","security_group":"sg-0test",
 "subnet_ids":["subnet-a","subnet-b"],"repo_url":"https://github.com/aburan28/crypto.git",
 "default_ttl":"4h","max_ttl_hours":48,"managed_tag":"ephemeral:managed",
 "cpu_ami_param":"/cpu","arm_ami_param":"/arm","gpu_ami_param":"/gpu","instance_profile":"ephemeral-runner"}
EOF

pass=0 fail=0
check() { # description, then a command that must succeed
    local what=$1
    shift
    if "$@"; then
        echo "  ok   $what"
        pass=$((pass + 1))
    else
        echo "  FAIL $what"
        fail=$((fail + 1))
    fi
}
last() { cat "$XDG_CACHE_HOME/ephemeral/last"; }
job_dir() { echo "$FAKE_AWS/s3/test-runs/jobs/$1"; }
state_of() { # job -> fake instance state
    local f j st
    for f in "$FAKE_AWS"/ec2/i-*; do
        case ${f##*/} in *.*) continue ;; esac
        read -r j st <"$f"
        [ "$j" = "$1" ] && echo "$st"
    done
    return 0
}
wait_state() { # job, instance state, seconds
    local i
    for ((i = 0; i < $3; i++)); do
        [ "$(state_of "$1")" = "$2" ] && return 0
        sleep 1
    done
    return 1
}
wait_phase() { # job, phase, seconds
    local i
    for ((i = 0; i < $3; i++)); do
        grep -q "\"phase\":\"$2\"" "$(job_dir "$1")/status.json" 2>/dev/null && return 0
        sleep 1
    done
    return 1
}

cd "$ROOT"

echo "== argument checks"
check "--ttl over the ceiling is refused" bash -c "! '$CLI' run --ttl 3d --no-source -- true 2>/dev/null"
check "--ttl under 15m is refused" bash -c "! '$CLI' run --ttl 5m --no-source -- true 2>/dev/null"
check "a missing command is refused" bash -c "! '$CLI' run --no-source 2>/dev/null"
check "--ref with --no-source is refused" bash -c "! '$CLI' run --ref HEAD --no-source -- true 2>/dev/null"
check "an unpushed --ref is refused" bash -c "
    c=\$(git commit-tree -m unpushed HEAD^{tree} -p HEAD)
    ! '$CLI' run --ref \$c -- true 2>/dev/null"
check "no instance was launched by the checks above" bash -c "! grep -q 'ec2 run-instances' '$FAKE_AWS/calls.log' 2>/dev/null"

echo "== working tree + sparse checkout + env + artifacts (exit 0)"
set +e
"$CLI" run --toolchain none --ttl 30m --paths taskq --env 'GREETING=hello world' \
    --artifact taskq/pyproject.toml --name e2e-ok -- \
    'echo "greeting=$GREETING"; test -f taskq/deploy/ephemeral/runner.sh && echo patch-applied; nproc >"$EPHEMERAL_OUT/nproc.txt"; ls'
rc=$?
set -e
id=$(last)
out=$XDG_CACHE_HOME/ephemeral/runs/$id
check "run exits 0" test "$rc" -eq 0
check "status is done" grep -q '"phase":"done"' "$out/status.json"
check "env var reached the job" grep -q 'greeting=hello world' "$out/job.log"
check "local (untracked) changes were applied" grep -q patch-applied "$out/job.log"
check "sparse checkout fetched only taskq" bash -c "! grep -qx src '$out/job.log' && grep -qx taskq '$out/job.log'"
check "EPHEMERAL_OUT came back" test -s "$out/out/nproc.txt"
check "--artifact came back" test -s "$out/out/artifacts/taskq/pyproject.toml"
check "host record and rusage uploaded" test -s "$out/host.txt" -a -s "$out/rusage.txt"
check "deadline timer was armed first" grep -q 'on-active=' "$FAKE_AWS"/ec2/*.timer
check "instance terminated itself" wait_state "$id" terminated 30
check "shutdown was called" bash -c "cat '$FAKE_AWS'/ec2/*.shutdown | grep -q 'shutdown -h now'"
check "deadline tag set at launch" grep -q "ephemeral:deadline,Value=" "$FAKE_AWS/calls.log"

echo "== failing command, argv form, capacity retry"
rm -f "$FAKE_AWS/capacity-fails"
set +e
FAKE_CAPACITY_FAIL=1 "$CLI" run --toolchain none --ttl 30m --no-source --name e2e-fail -- \
    bash -c 'printf "%s\n" "two words"; exit 3'
rc=$?
set -e
id=$(last)
out=$XDG_CACHE_HOME/ephemeral/runs/$id
check "run exits with the command's status (3)" test "$rc" -eq 3
check "status is failed with exit_code 3" bash -c "grep -q '\"phase\":\"failed\"' '$out/status.json' && grep -q '\"exit_code\":3' '$out/status.json'"
check "argv quoting kept 'two words' together" grep -qx 'two words' "$out/job.log"
check "retried the next subnet after a capacity error" test "$(grep -c 'ec2 run-instances' "$FAKE_AWS/calls.log")" -ge 3
check "failed job still terminated" wait_state "$id" terminated 30
check "fetched log is complete (final upload before final status)" grep -q terminating "$out/job.log"

echo "== detach, then kill"
id=$("$CLI" run --toolchain none --ttl 30m --no-source --detach --name e2e-kill -- sleep 600)
check "--detach printed the job id" test -n "$id"
check "job reached running" wait_phase "$id" running 180
"$CLI" kill "$id" >/dev/null 2>&1
check "kill marked the instance terminated" test "$(state_of "$id")" = terminated
sleep 10
check "status stays killed after the runner exits" grep -q '"phase":"killed"' "$(job_dir "$id")/status.json"
check "ls lists the jobs" bash -c "'$CLI' ls | grep -q e2e-kill"
check "logs prints the log" bash -c "'$CLI' logs '$id' | grep -q 'phase: running'"

echo
echo "$pass passed, $fail failed"
[ "$fail" -eq 0 ]
