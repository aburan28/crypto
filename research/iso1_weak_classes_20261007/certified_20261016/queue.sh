#!/bin/bash
cd /home/user/iso1_cert
for p in 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199; do
  [ -s p${p}_cert.csv.gz ] && continue
  echo -e "$(date +%s)\t$p\tRUNNING" >> status.tsv
  if RAYON_NUM_THREADS=4 ./iso1_class_census_cert $p p${p}_cert.csv --derive-twists --absolute-orbit-quotient --certified 2> p${p}_cert.stderr; then
    sha256sum p${p}_cert.csv > p${p}_cert.csv.sha256
    gzip -9 p${p}_cert.csv
    echo -e "$(date +%s)\t$p\tCOMPLETE" >> status.tsv
  else
    echo -e "$(date +%s)\t$p\tFAILED exit=$?" >> status.tsv; exit 1
  fi
done
