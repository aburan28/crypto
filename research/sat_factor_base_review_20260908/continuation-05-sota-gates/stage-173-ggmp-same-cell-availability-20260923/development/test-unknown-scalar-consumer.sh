#!/bin/zsh
set -eu

repo=/Volumes/SSD990/crypto-f4-single-target-stage171
backup=/Volumes/SSD990/kic-stage173-ggmp-dev/Cargo.lock.test-backup
/bin/cp "$repo/Cargo.lock" "$backup"
restore_lock() {
  /bin/cp "$backup" "$repo/Cargo.lock"
}
trap restore_lock EXIT
/bin/cp "$repo/research/weil_factor_composition_20260914/validation/dependencies.lock.txt" "$repo/Cargo.lock"
cd "$repo"
python3 scripts/test_koblitz_unknown_scalar_panel.py -q
