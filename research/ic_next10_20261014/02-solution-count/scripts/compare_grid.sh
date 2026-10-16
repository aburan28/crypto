#!/usr/bin/env bash
# Compare the Rust group-law census (per-cell histogram) with the gp direct grid count (per target).
# Usage: compare_grid.sh <validation dir>; prints arm,m,t,rust_hist,gp_hist,match
V=$1
awk -F, 'NR==FNR { if (FNR>1) { key=$1","$2","$3; c[key]=c[key]" "$6 } ; next }
  FNR>1 { key=$5","$6","$7; if (!(key in c)) next;
    n=split(c[key], a, " "); for(i=2;i<=n;i++){x=a[i]+0;j=i-1;while(j>=1 && a[j]+0>x){a[j+1]=a[j];j--};a[j+1]=x}
    s=""; i=1; while(i<=n){j=i; while(j<n && a[j+1]+0==a[i]+0) j++; s=s (s==""?"":" ") a[i]":"(j-i+1); i=j+1}
    print key","$20","s","($20==s?"MATCH":"DIFF") }' $V/grid_counts.csv $V/rust_counts.csv
