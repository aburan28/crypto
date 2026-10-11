#!/bin/bash
set -euo pipefail
export PATH=/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin
cd /Volumes/SSD990/crypto/worktrees/iso1-theorem-20261009
exec research/iso1_weak_classes_20261007/census_runtime/target/release/iso1_census_queue \
  '/Users/adamburan/Library/Application Support/crypto-iso1/census-20261009-run1' \
  59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131 137 139 149 151 157 163 167 173 179 181 191 193 197 199
