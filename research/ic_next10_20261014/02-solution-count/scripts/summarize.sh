#!/usr/bin/env bash
# Derive per-cell ratios from data/census.csv (awk only; no new computation of counts).
# Output data/ratios.csv: one row per (curve, arm, m, t) with
#   r_rand = mean / mean_random (same curve, m, t, targets), r_exact = mean / ((2 B_E)^m / n),
#   r_product = mean / (B^m / p), r_Bm1 = mean / B^(m-1), disp = var/mean, disp_rand.
cd "$(dirname "$0")/.."
awk -F, 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
 { key=$h["p"]","$h["m"]","$h["t"]; if($h["arm"]=="random"){rm[key]=$h["mean"]; rv[key]=$h["var"]} row[NR]=$0 }
 END{ print "bits,sign,p,arm,m,t,B,B_E,targets,mean,ci_lo,ci_hi,mean_random,r_rand,r_exact,r_product,r_Bm1,disp,disp_rand,frac_ge1,scored";
  for(r=2;r<=NR;r++){ split(row[r],f,","); key=f[h["p"]]","f[h["m"]]","f[h["t"]];
   mr=rm[key]; mean=f[h["mean"]];
   rr=(mr>0)?mean/mr:"NA"; re=(f[h["pred_exact"]]>0)?mean/f[h["pred_exact"]]:"NA";
   rp=mean/f[h["pred_product"]]; rb=mean/f[h["Bm1"]];
   d=(mean>0)?f[h["var"]]/mean:"NA"; dr=(mr>0)?rv[key]/mr:"NA";
   sc=(mr>=0.5)?1:0;
   printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d\n", f[1],f[2],f[h["p"]],f[h["arm"]],f[h["m"]],f[h["t"]],f[h["B"]],f[h["B_E"]],f[h["targets"]],mean,f[h["ci_lo"]],f[h["ci_hi"]],mr,rr,re,rp,rb,d,dr,f[h["frac_ge1"]],sc } }' data/census.csv > data/ratios.csv
