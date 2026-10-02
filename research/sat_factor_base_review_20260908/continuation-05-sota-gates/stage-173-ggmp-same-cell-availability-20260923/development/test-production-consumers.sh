#!/bin/zsh
set -eu

repo=/Volumes/SSD990/crypto-f4-single-target-stage171
cd "$repo"
python3 scripts/test_koblitz_stage23_terminal_evidence.py -q
python3 research/sat_factor_base_review_20260908/continuation-05-sota-gates/verify_stage18_degree23_panel.py --self-test
