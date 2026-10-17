#!/bin/sh
# External reference (AGENTS.md §15): PARI/GP `elllog` on the study's curves
# and targets, single-threaded, timed by gp's own getabstime() around the
# elllog call only.  Thin orchestration: the shell substitutes the frozen
# inputs into a gp program; gp checks #E = n before timing, and prints the
# recovered log for the caller to check against the native arm.
#
#   pari_elllog.sh <raw_dir> <role> <bits> <k> <label>
#
# <role> is bsgs or rho (the curve's role in inputs/curves.json); <k> is the
# known answer the native arm used (from its JSON line: rho prints `k`, and the
# BSGS driver's one-target `answers_digest` inverts to it, see the report).
set -eu
raw=$1; role=$2; bits=$3; k=$4; label=$5
here=$(cd "$(dirname "$0")" && pwd)
curves="$here/../inputs/curves.json"
c=$(jq -c ".curves[] | select(.role == \"$role\" and .bits == $bits)" "$curves")
p=$(echo "$c" | jq -r .p); a=$(echo "$c" | jq -r .a); b=$(echo "$c" | jq -r .b)
gx=$(echo "$c" | jq -r .generator.x); gy=$(echo "$c" | jq -r .generator.y)
n=$(echo "$c" | jq -r .subgroup_order)
gp -q -f -s 512M <<EOF >> "$raw/pari.jsonl"
default(nbthreads, 1);
E = ellinit([$a, $b], $p);
if (ellcard(E) != $n, error("group order mismatch"));
G = [Mod($gx, $p), Mod($gy, $p)];
Q = ellmul(E, G, $k);
t0 = getabstime(); x = elllog(E, Q, G, $n); t = getabstime() - t0;
printf("{\"label\":\"%s\",\"tool\":\"PARI/GP %s\",\"role\":\"%s\",\"bits\":%d,\"k\":%d,\"log\":%d,\"verified\":%s,\"elllog_ms\":%d}\n", "$label", strjoin(apply(x->Str(x), Vec(version())[1..3]), "."), "$role", $bits, $k, x, if (x == $k && ellmul(E, G, x) == Q, "true", "false"), t);
EOF
