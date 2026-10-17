#!/bin/bash
# interleaved: blocked (B=1), control (B=0), msolve; 1 thread; 3 rounds
cd "$(dirname "$0")"
out=ab_full.tsv; : > $out
cases=$(GB_CAS_LARGE=1 ./gce_blk fp inL 1 2>/dev/null | sed -n 's/.*"case":"\([^"]*\)".*/\1/p')
for r in 1 2 3; do
  for c in $cases; do
    for b in 1 0; do
      F4_FP_FL_BLOCKED=$b GB_CAS_LARGE=1 GB_CAS_ONLY=$c RAYON_NUM_THREADS=1 ./gce_blk fp inL 1 \
        | python3 -c "import sys,json; d=json.loads(sys.stdin.read()); print('$r\t$c\tours_B$b\t%.3f\t%s' % (d['ours_f4_ms_median'], d['basis_fp']))" >> $out
    done
    a=$(date +%s%N); msolve -g 2 -t 1 -f inL/$c.ms -o /dev/null >/dev/null 2>&1; b=$(date +%s%N)
    printf '%s\t%s\tmsolve\t%s\t-\n' $r $c $(( (b-a)/1000 ))e-3 >> $out
  done
done
echo done
