#!/bin/sh
# lit2026 W11 (b): interleaved end-to-end Koblitz class walk, (n, a2) = (19, 1)
# (ℓ = 457 Vélu steps), original / previous / candidate (thin orchestration).
#
#   run_classwalk_arms.sh <outdir> <rounds> <orig_bin> <cand_bin> <repo_root>
#
# Each run happens in its own scratch cwd holding a copy of
# experiments/koblitz_isogeny_cost_sweep.json (the census the walk checks
# against), so the committed snapshot is never overwritten.  The snapshot
# each run writes is kept, and the arms' snapshots are compared with their
# "seconds" fields removed.
set -eu
outdir=$1; rounds=$2; orig=$3; cand=$4; repo=$5
mkdir -p "$outdir"
: > "$outdir/samples.tsv"
arms="original previous candidate"
r=1
while [ "$r" -le "$rounds" ]; do
  set -- $arms
  i=0
  while [ $i -lt $(( (r - 1) % 3 )) ]; do
    first=$1; shift; set -- "$@" "$first"; i=$((i + 1))
  done
  for arm in "$@"; do
    wd="$outdir/wd-$arm-$r"
    mkdir -p "$wd/experiments"
    cp "$repo/experiments/koblitz_isogeny_cost_sweep.json" "$wd/experiments/"
    up=$(uptime | sed -E 's/.*load averages?: //')
    t0=$(perl -MTime::HiRes=time -e 'printf "%.6f", time')
    case $arm in
      original) (cd "$wd" && "$orig" 19:1 > stdout.txt 2>&1) ;;
      previous) (cd "$wd" && env -u BTW_XEVAL "$cand" 19:1 > stdout.txt 2>&1) ;;
      candidate) (cd "$wd" && BTW_XEVAL=fraction "$cand" 19:1 > stdout.txt 2>&1) ;;
    esac
    t1=$(perl -MTime::HiRes=time -e 'printf "%.6f", time')
    printf '%s\t%s\t%s\t%s\n' "$arm" "$r" "$(perl -e "printf '%.3f', $t1 - $t0")" "$up" >> "$outdir/samples.tsv"
  done
  r=$((r + 1))
done
