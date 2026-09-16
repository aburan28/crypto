#!/usr/bin/env bash
# One matched set of EC2 spot quotes for the four instance types the campaign
# can choose between: g6 (L4), g6e (L40S), g7 (RTX PRO 4500), g7e (RTX PRO 6000).
#
# Matched is the point. An Ada quote taken on one day and a Blackwell quote
# taken on another are not a comparison -- spot moves, and ../../RTX-PRO4500.md
# records a g7-vs-g7e break-even only a few percent wide, narrow enough for a
# day's drift to flip it. So take all four at once, from one invocation, and
# record the timestamp with them.
#
#   bash benchmarks/ada/prices.sh                     # us-west-2, us-east-1, us-east-2
#   REGIONS="us-west-2" SIZE=4xlarge bash benchmarks/ada/prices.sh
#
# Writes benchmarks/ada/prices.json and prints a markdown table for
# ../../ADA-L4-L40S.md. Needs only the AWS CLI with describe-spot-price-history
# permission; it reads and never launches anything.
set -uo pipefail
cd "$(dirname "$0")/../.." || exit 1
OUT="${OUT:-benchmarks/ada/prices.json}"
REGIONS="${REGIONS:-us-west-2 us-east-1 us-east-2}"
SIZE="${SIZE:-2xlarge}"
TYPES="${TYPES:-g6.$SIZE g6e.$SIZE g7.$SIZE g7e.$SIZE}"

if ! command -v aws >/dev/null 2>&1; then
  echo "REFUSING: no aws CLI on PATH, so no price can be observed." >&2
  echo "Install it, or run this from a host that has it; do not fill the table" >&2
  echo "in ../../ADA-L4-L40S.md from anywhere else." >&2
  exit 2
fi

python3 - "$OUT" "$REGIONS" "$TYPES" <<'PY'
import json, subprocess, sys, time

out_path, regions, types = sys.argv[1], sys.argv[2].split(), sys.argv[3].split()
GPU = {"g6": "L4", "g6e": "L40S", "g7": "RTX PRO 4500", "g7e": "RTX PRO 6000"}

rec = {"takenAt": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
       "method": "describe-spot-price-history, most recent Linux/UNIX quote, "
                 "cheapest availability zone per type; one invocation so every "
                 "quote is from the same moment",
       "regions": regions, "types": types, "quotes": [], "errors": []}

for region in regions:
    for itype in types:
        cmd = ["aws", "ec2", "describe-spot-price-history",
               "--region", region, "--instance-types", itype,
               "--product-descriptions", "Linux/UNIX",
               "--start-time", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "--query", "SpotPriceHistory[].[AvailabilityZone,SpotPrice,Timestamp]",
               "--output", "json"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            # A type not offered in a region is an ordinary answer, not a
            # failure of the run. Record it and keep going; a missing row is
            # never a zero price.
            rec["errors"].append({"region": region, "instanceType": itype,
                                  "returncode": p.returncode,
                                  "stderr": p.stderr.strip()[-600:]})
            continue
        rows = json.loads(p.stdout or "[]")
        if not rows:
            rec["errors"].append({"region": region, "instanceType": itype,
                                  "reason": "no spot price history returned"})
            continue
        az, price, stamp = min(rows, key=lambda r: float(r[1]))
        rec["quotes"].append({"region": region, "instanceType": itype,
                              "gpu": GPU.get(itype.split(".")[0], "?"),
                              "availabilityZone": az,
                              "usdPerHour": float(price), "quotedAt": stamp})

json.dump(rec, open(out_path, "w"), indent=1)
print(f"wrote {out_path}\n")
print("| region | type | GPU | $/hr |")
print("|---|---|---|---:|")
for q in sorted(rec["quotes"], key=lambda q: (q["region"], q["instanceType"])):
    print("| %s | %s | %s | %.4f |" % (q["region"], q["instanceType"],
                                       q["gpu"], q["usdPerHour"]))
if rec["errors"]:
    print("\nno quote (type not offered there, or the call failed):")
    for e in rec["errors"]:
        print("  %s %s: %s" % (e["region"], e["instanceType"],
                               e.get("reason") or e.get("stderr", "")))
if not rec["quotes"]:
    sys.exit(1)      # an empty table is not a price observation
PY
