#!/usr/bin/env bash
set -euo pipefail
cd /workspace/repos/crypto
python3 scripts/build_lab_browser.py
PYTHONPATH=scripts/site python3 -m unittest test_build.BuildTests.test_lab_browser_index_is_current_and_self_consistent
python3 scripts/test_check_ic_leaderboard_paths.py
python3 scripts/check_ic_leaderboard_paths.py
python3 scripts/build_ic_leaderboard.py --check
node tools/render_ic_dashboard.mjs --check
python3 -m py_compile research/full_volcano_graphwide_ic_20261017/execution/render_survey_checkpoint.py research/full_volcano_graphwide_ic_20261017/execution/render_original_local_strata.py research/full_volcano_graphwide_ic_20261017/execution/package_survey_checkpoint.py
mkdir /workspace/volcano_recovered_viewer_2935_20261017
/workspace/scratch/volcano-graph-export /workspace/complete_volcano_graphs_20261017 /workspace/volcano_ic_repaired_checkpoint_20261017T1450Z/inventory.json research/full_volcano_graphwide_ic_20261017/graph-viewer/template.html /workspace/volcano_recovered_viewer_2935_20261017/graphs.html > /workspace/volcano_recovered_viewer_2935_20261017/export.json
node research/full_volcano_graphwide_ic_20261017/graph-viewer/test-ui.mjs /workspace/volcano_recovered_viewer_2935_20261017/graphs.html > /workspace/volcano_recovered_viewer_2935_20261017/ui.log
python3 -m http.server 8774 --bind 127.0.0.1 --directory /workspace/volcano_recovered_viewer_2935_20261017 > /workspace/volcano_recovered_viewer_2935_20261017/http.log 2>&1 &
viewer_http_pid=$!
trap 'kill "$viewer_http_pid"' EXIT
sleep 1
python3 research/full_volcano_graphwide_ic_20261017/graph-viewer/capture.py http://127.0.0.1:8774/graphs.html /workspace/volcano_recovered_viewer_2935_20261017/screenshots > /workspace/volcano_recovered_viewer_2935_20261017/chromium.log 2>&1
