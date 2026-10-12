#!/usr/bin/env bash
# one.sh SUITE RUNG ARM OUTDIR : callgrind one rung, collection limited to groebner_decompose
set -uo pipefail
suite=$1; rung=$2; arm=$3; out=$4
bench=/home/user/crypto/target/release/examples/groebner_stage_bench
ladder=(); [ "$suite" != frozen ] && ladder=(--ladder "$suite")
d="$out/$suite/$arm/rung$rung"; mkdir -p "$d"
[ -e "$d/summary.json" ] && exit 0
KIC_F4_MULTIPLIERS=$arm KIC_CHAIN_ORDER=interleaved KIC_LINEAR_ELIM=1 KIC_F4_DROP=complete \
  valgrind --tool=callgrind --toggle-collect='*groebner_decompose*' --callgrind-out-file="$d/callgrind.out" \
  "$bench" --label "$arm" "${ladder[@]}" --rung "$rung" --out "$d/stage" > "$d/stdout" 2> "$d/valgrind.log"
if ! python3 -c "import json,sys; sys.exit(0 if json.load(open('$d/stage/stage.json'))['rows'] else 1)" 2>/dev/null; then
  echo '{"skipped": true}' > "$d/summary.json"; rm -f "$d/callgrind.out"; exit 0
fi
callgrind_annotate --inclusive=yes --threshold=100 "$d/callgrind.out" 2>/dev/null > "$d/inclusive.txt"
python3 - "$d" <<'PY'
import json, re, sys, pathlib
d = pathlib.Path(sys.argv[1])
text = (d / "inclusive.txt").read_text().splitlines()
total = None; funcs = {}
for line in text:
    m = re.match(r"\s*([\d,]+)\s+\(\s*[\d.]+%\)\s+(PROGRAM TOTALS|\S.*)$", line)
    if not m: continue
    ir = int(m.group(1).replace(",", "")); name = m.group(2)
    if name == "PROGRAM TOTALS": total = ir; continue
    name = re.sub(r"\s*\[.*\]$", "", name).replace("???:", "")
    funcs[name] = max(funcs.get(name, 0), ir)
row = json.loads((d / "stage" / "stage.json").read_text())["rows"][0]
keep = {k: v for k, v in funcs.items() if v >= 0.002 * total}
json.dump({"rung": row["curve"] + f" m={row['m']} fi={row['factor_index']} first={row.get('first_target',0)}",
           "word_ops": row["word_ops"], "specialise_word_ops": row["specialise_word_ops"],
           "reductions": row["reductions"], "ir_total": total, "inclusive_ir": keep},
          open(d / "summary.json", "w"), indent=1)
PY
gzip -9 -f "$d/callgrind.out"
