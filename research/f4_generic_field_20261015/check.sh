#!/usr/bin/env bash
# Gates G2 and G4 (PROTOCOL.md §3) on the committed rows; writes runs/checks.txt.
# Usage: check.sh CHECKER_BIN
set -uo pipefail
C=$1
cd "$(dirname "$0")/runs"
tmp=$(mktemp -d)
{
  echo "# G2: u32 engine vs generic Fp32 (engine label rewritten; backend, builder ignored)"
  for s in kummer-p2 dickson-p1 torus-T31 torus-T24; do
    sed 's/"engine":"f4_tower_generic"/"engine":"f4_fp_tower"/' "eq/$s.fp32.jsonl" > "$tmp/$s.jsonl"
    echo "## $s: $("$C" compare-rows "eq/$s.u32.jsonl" "$tmp/$s.jsonl" --ignore backend,builder 2>&1 | tail -1)"
  done
  echo "# G2: Kummer p2 at Fp32 vs Fp64, Fp52, FpWide<4> (backend ignored)"
  for b in fp64 fp52 wide4; do
    echo "## $b: $("$C" compare-rows eq/kummer-p2.fp32.jsonl "eq/kummer-p2.$b.jsonl" --ignore backend 2>&1 | tail -1)"
  done
  echo "# G4: pkm_tower_check verify on the u64-builder depth rows"
  for f in depth/*.jsonl ext/*.jsonl; do
    if grep -q '"builder":"u64"' "$f"; then
      echo "## $f"
      "$C" verify "$f" 2>&1 | tail -3
    fi
  done
  echo "# G4: wide-builder inline verification (finished tower rows)"
  for f in depth/*.jsonl ext/*.jsonl; do
    if grep -q '"builder":"wide"' "$f"; then
      n=$(grep -c '"verify":{' "$f" || true)
      ok=$(grep -c '"ok":true' "$f" || true)
      bad=$(grep -c '"ok":false' "$f" || true)
      echo "## $f: verified rows $n, ok $ok, failed $bad"
    fi
  done
  echo "# planted_ok false anywhere:"
  grep -l '"planted_ok":false' depth/*.jsonl eq/*.jsonl || echo "none"
} > checks.txt
rm -rf "$tmp"
cat checks.txt
