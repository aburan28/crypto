#!/bin/sh
set -eu
cd "$(dirname "$0")"
tmp="published_cases.sqlite.tmp.$$"
trap 'rm -f "$tmp"' EXIT HUP INT TERM
sqlite3 -bail "$tmp" < schema.sql
sqlite3 -bail "$tmp" <<'SQL'
PRAGMA foreign_keys = ON;
.mode csv
.import --skip 1 cases.csv published_case
.import --skip 1 minimum_case_patterns.csv minimum_case_pattern
CREATE TEMP TABLE assertion (pass INTEGER NOT NULL CHECK (pass = 1));
INSERT INTO assertion SELECT CASE WHEN count(*) = 36 AND min(case_no) = 1 AND max(case_no) = 36
  AND sum(CASE WHEN e = 0 THEN 1 ELSE 0 END) = 5
  AND sum(CASE WHEN e > 0 THEN 1 ELSE 0 END) = 31
  THEN 1 ELSE 0 END
FROM published_case;
INSERT INTO assertion SELECT CASE WHEN count(*) = 5 AND
  count(*) = (SELECT count(*) FROM published_case WHERE e = 0)
  THEN 1 ELSE 0 END FROM minimum_case_recipe;
SELECT 'case-index-ok';
PRAGMA integrity_check;
SQL
mv "$tmp" published_cases.sqlite
trap - EXIT HUP INT TERM
