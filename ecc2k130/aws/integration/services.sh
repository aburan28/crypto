#!/usr/bin/env bash
# Bring up Postgres, Redis and an S3 endpoint for test_integration.py.
#
#   services.sh up      start them, wait until they answer, write the env file
#   services.sh down    stop them and take their data with them
#   services.sh env     print the env file's path (or the exports themselves)
#
# Docker is the intended way and the one CI uses.  When there is no daemon --
# a locked-down workstation, a container without a socket -- the same services
# are started from the binaries on the host instead, because a suite that only
# runs in one place is a suite that rots in the other.  Either way the tests
# see nothing but three URLs.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"           # ecc2k130/
state="$root/build/integration"
envfile="$state/env"

PG_PORT="${ECC_IT_PG_PORT:-5433}"
REDIS_PORT="${ECC_IT_REDIS_PORT:-6380}"
S3_PORT="${ECC_IT_S3_PORT:-9010}"
BUCKET="${ECC_IT_BUCKET:-ecc2k130-integration}"

say() { printf '%s\n' "$*" >&2; }
die() { say "services.sh: $*"; exit 1; }

haveDocker() {
    [ "${ECC_IT_NATIVE:-0}" = "1" ] && return 1
    command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1
}

compose() { docker compose -f "$here/docker-compose.yml" "$@"; }

# Wait for a TCP port to answer at all.  Every service gets a protocol-level
# check after this one; this is only so the protocol check is not racing the
# container's first second.
waitPort() {
    local host=$1 port=$2 name=$3 tries=${4:-150}
    for _ in $(seq "$tries"); do
        if (exec 3<>"/dev/tcp/$host/$port") 2>/dev/null; then exec 3<&- 3>&-; return 0; fi
        sleep 0.2
    done
    die "$name did not open $host:$port"
}

portOpen() {
    (exec 3<>"/dev/tcp/${1}/${2}") 2>/dev/null && { exec 3<&- 3>&-; return 0; }
    return 1
}

waitHttp() {
    local url=$1 name=$2 tries=${3:-150}
    for _ in $(seq "$tries"); do
        if curl -fsS -o /dev/null "$url" 2>/dev/null; then return 0; fi
        sleep 0.2
    done
    die "$name did not answer $url"
}

checkClients() {
    python3 - <<'PY' || die "missing python clients: python3 -m pip install -r $here/requirements.txt"
import importlib.util as util
import sys
missing = [m for m in ("psycopg", "redis", "boto3") if util.find_spec(m) is None]
if missing:
    print("missing: %s" % ", ".join(missing), file=sys.stderr)
sys.exit(1 if missing else 0)
PY
}

writeEnv() {
    local url=$1
    mkdir -p "$state"
    cat > "$envfile" <<ENVEOF
# Written by services.sh; read by \`make test-integration\`.
ECC_INTEGRATION=1
RHO_CAMPAIGN=integration
DATABASE_URL=$url
RHO_REDIS_URL=redis://127.0.0.1:$REDIS_PORT/0
ECC_BUCKET=$BUCKET
ECC_S3_ENDPOINT=http://127.0.0.1:$S3_PORT
ECC_S3_REGION=us-east-1
AWS_ACCESS_KEY_ID=$2
AWS_SECRET_ACCESS_KEY=$3
AWS_DEFAULT_REGION=us-east-1
ENVEOF
    say "services.sh: environment in $envfile"
}

# -- docker ----------------------------------------------------------------
upDocker() {
    say "services.sh: starting postgres, redis and minio under docker"
    ECC_IT_PG_PORT="$PG_PORT" ECC_IT_REDIS_PORT="$REDIS_PORT" ECC_IT_S3_PORT="$S3_PORT" \
        compose up -d --wait --wait-timeout 120 postgres redis
    ECC_IT_PG_PORT="$PG_PORT" ECC_IT_REDIS_PORT="$REDIS_PORT" ECC_IT_S3_PORT="$S3_PORT" \
        compose up -d minio
    waitPort 127.0.0.1 "$S3_PORT" minio
    waitHttp "http://127.0.0.1:$S3_PORT/minio/health/live" minio
    writeEnv "postgresql://controlplane:integration@127.0.0.1:$PG_PORT/controlplane" \
             ecc2k130 ecc2k130-integration
}

downDocker() {
    ECC_IT_PG_PORT="$PG_PORT" ECC_IT_REDIS_PORT="$REDIS_PORT" ECC_IT_S3_PORT="$S3_PORT" \
        compose down -v --remove-orphans || true
}

# -- native ----------------------------------------------------------------
pgBin() {
    local exe=$1
    if command -v "$exe" >/dev/null 2>&1; then command -v "$exe"; return; fi
    local found
    found="$(ls -d /usr/lib/postgresql/*/bin 2>/dev/null | sort -V | tail -1 || true)"
    [ -n "$found" ] && [ -x "$found/$exe" ] && { printf '%s\n' "$found/$exe"; return; }
    die "no $exe on PATH and no /usr/lib/postgresql/*/bin (install postgresql, or start docker)"
}

# initdb and postgres refuse to run as root, so when this is root the data
# directory is handed to the postgres account the package created.
asPostgres() {
    if [ "$(id -u)" = "0" ] && id postgres >/dev/null 2>&1; then
        su postgres -c "$1"
    else
        bash -c "$1"
    fi
}

upNative() {
    say "services.sh: no docker daemon; starting postgres, redis and moto from the host"
    command -v redis-server >/dev/null 2>&1 || die "no redis-server (install redis, or start docker)"
    command -v moto_server  >/dev/null 2>&1 || \
        die "no moto_server for the S3 endpoint: python3 -m pip install 'moto[server]'"
    local initdb pgctl data
    initdb="$(pgBin initdb)"; pgctl="$(pgBin pg_ctl)"
    data="$state/pg"
    mkdir -p "$state/redis" "$state/run"
    if portOpen 127.0.0.1 "$PG_PORT"; then
        say "services.sh: postgres already answering on $PG_PORT; leaving it"
    else
        if [ ! -s "$data/PG_VERSION" ]; then
            rm -rf "$data"; mkdir -p "$data"
            [ "$(id -u)" = "0" ] && id postgres >/dev/null 2>&1 && chown postgres "$state" "$data"
            asPostgres "$initdb -D '$data' -U postgres --auth=trust -E UTF8 >/dev/null"
        fi
        asPostgres "$pgctl -D '$data' -o '-p $PG_PORT -k \"$data\" -c listen_addresses=127.0.0.1 -c fsync=off -c max_connections=200' -l '$state/postgres.log' -w start" \
            || { say "--- postgres.log"; tail -20 "$state/postgres.log" >&2 || true; die "postgres did not start"; }
        waitPort 127.0.0.1 "$PG_PORT" postgres
    fi

    if portOpen 127.0.0.1 "$REDIS_PORT"; then
        say "services.sh: redis already answering on $REDIS_PORT; leaving it"
    else
        redis-server --port "$REDIS_PORT" --bind 127.0.0.1 --save '' --appendonly no \
            --daemonize yes --dir "$state/redis" --pidfile "$state/run/redis.pid" \
            --logfile "$state/redis.log"
        waitPort 127.0.0.1 "$REDIS_PORT" redis
    fi

    if portOpen 127.0.0.1 "$S3_PORT"; then
        say "services.sh: an S3 endpoint already answering on $S3_PORT; leaving it"
    else
        moto_server -p "$S3_PORT" -H 127.0.0.1 > "$state/moto.log" 2>&1 &
        echo $! > "$state/run/moto.pid"
        waitPort 127.0.0.1 "$S3_PORT" "the S3 endpoint"
    fi
    waitHttp "http://127.0.0.1:$S3_PORT/" "the S3 endpoint"

    writeEnv "postgresql://postgres@127.0.0.1:$PG_PORT/postgres" integration integration
}

downNative() {
    local pgctl data
    data="$state/pg"
    if [ -s "$data/PG_VERSION" ]; then
        pgctl="$(pgBin pg_ctl)"
        asPostgres "$pgctl -D '$data' -m immediate -w stop" >/dev/null 2>&1 || true
    fi
    [ -s "$state/run/redis.pid" ] && kill "$(cat "$state/run/redis.pid")" 2>/dev/null || true
    [ -s "$state/run/moto.pid" ]  && kill "$(cat "$state/run/moto.pid")"  2>/dev/null || true
    rm -rf "$state/pg" "$state/redis" "$state/run"
}

case "${1:-up}" in
    up)
        checkClients
        if haveDocker; then upDocker; else upNative; fi
        ;;
    down)
        if haveDocker; then downDocker; else downNative; fi
        rm -f "$envfile"
        ;;
    env)
        [ -s "$envfile" ] || die "no environment yet: services.sh up"
        if [ "${2:-}" = "--export" ]; then sed 's/^\([A-Z]\)/export \1/' "$envfile"
        else printf '%s\n' "$envfile"; fi
        ;;
    *)
        die "usage: services.sh [up|down|env]"
        ;;
esac
