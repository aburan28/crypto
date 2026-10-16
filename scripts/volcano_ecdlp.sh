#!/usr/bin/env bash
# volcano_ecdlp.sh - point at a curve, measure ECDLP across its isogeny volcano.
#
#   scripts/volcano_ecdlp.sh --curve P,A,B [--order N] [volcano-rho options]
#   scripts/volcano_ecdlp.sh --ca [PATH] --curve P,A,B ...   rho by cryptanalysis's `ca`
#   scripts/volcano_ecdlp.sh --search [volcano-rho options]  let it find a curve
#
# Builds `volcano-rho` (release), then runs it on y^2 = x^3 + Ax + B over F_P
# (P prime < 2^62). It counts points unless --order is given, draws a seeded
# random generator and secret k (--seed), walks the ell-isogeny volcano
# (--ell 2|3) with verified Vélu maps, carries the same (P, kP) to every
# curve, and runs Pollard rho --trials times on each. It prints min / mean /
# max per curve and per level, and writes to OUT (default
# volcano-ecdlp-<UTC time>/):
#
#   results.json      curve, volcano, every summary
#   trials.csv        one row per solve
#   instances.jsonl   every curve's instance in `ca` format (p, a, b, order, g, h)
#   volcano.svg/.png, per_curve.*, per_level.*   when python3 has matplotlib
#
# Everything is deterministic code: the same command and --seed give the same
# curves, instance and iteration counts (wall times vary).
#
# --ca [PATH] uses `ca solve --alg rho` from aburan28/cryptanalysis as the rho
# engine (PATH defaults to $CA, then `ca` on PATH). Extra `ca` options go
# through --ca-arg, e.g. --ca-arg --no-negation.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
args=() out="" search=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ca)
            if [[ $# -gt 1 && "$2" != --* ]]; then ca="$2"; shift; else ca="${CA:-ca}"; fi
            args+=(--engine ca --ca "$ca") ;;
        --search) search=1 ;;
        --out) out="$2"; shift ;;
        -h|--help) sed -n '2,28p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) args+=("$1") ;;
    esac
    shift
done
if [[ $search -eq 0 ]] && ! printf '%s\n' "${args[@]+"${args[@]}"}" | grep -qx -- --curve; then
    echo "volcano_ecdlp.sh: pass --curve P,A,B (or --search to let it find one); --help for more" >&2
    exit 2
fi
out="${out:-volcano-ecdlp-$(date -u +%Y%m%dT%H%MZ)}"

cargo build --release --quiet --manifest-path "$root/Cargo.toml" --bin volcano-rho
"$root/target/release/volcano-rho" "${args[@]+"${args[@]}"}" --out "$out"

if python3 -c "import matplotlib" 2>/dev/null; then
    python3 "$root/research/isogeny-volcano-rho/plot.py" "$out" >/dev/null && echo "graphs: $out/volcano.svg, per_curve.svg, per_level.svg"
else
    echo "graphs skipped (python3 has no matplotlib); render later with research/isogeny-volcano-rho/plot.py $out"
fi
