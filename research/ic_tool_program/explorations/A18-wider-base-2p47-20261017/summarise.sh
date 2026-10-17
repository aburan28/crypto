#!/bin/bash
# summarise.sh pairs.tsv: per curve, the geometric mean of off/arm cold time
# over the clean pairs (neither run contended) and over all pairs, with a
# 95% t interval for the clean ones.
awk -F'\t' '
function tcrit(df) { split("12.706 4.303 3.182 2.776 2.571 2.447 2.365 2.306 2.262 2.228 2.201 2.179", t, " "); return t[df] }
NR > 1 { l = log($6); all[$1] += l; na[$1]++
         if ($7 == "false" && $8 == "false") { s[$1] += l; ss[$1] += l*l; n[$1]++ } }
END { printf "curve\tclean_pairs\toff_over_arm\tlo95\thi95\tall_pairs\tall_geomean\n"
      for (c in n) { m = s[c]/n[c]; v = (ss[c] - n[c]*m*m)/(n[c]-1); h = tcrit(n[c]-1)*sqrt(v/n[c])
        printf "%s\t%d\t%.3f\t%.3f\t%.3f\t%d\t%.3f\n", c, n[c], exp(m), exp(m-h), exp(m+h), na[c], exp(all[c]/na[c]) } }' "$1" | sort
