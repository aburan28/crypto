"""Custody and completeness audit for the signed P-192 short-support screen."""

import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
input_path=HERE/"type_i_line_guided_collector.json"
prior_path=HERE/"type_i_line_three_sum_coverage.json"
screen_path=HERE/"type_i_line_signed_short_relation_screen.json"
sage_path=HERE/"type_i_line_signed_short_sage_verification.json"
prior=json.loads(prior_path.read_text())
screen=json.loads(screen_path.read_text())
sage=json.loads(sage_path.read_text())
assert screen["status"]=="verified_p192_line_base_signed_short_relation_screen"
assert screen["source_sha256"]==sha(HERE/"line_signed_short_relation_screen.py")
assert screen["input_sha256"]==sha(input_path)
assert screen["prior_support_sha256"]==sha(prior_path)
assert screen["support_counts"]=={"1":141,"2":10011,"3":477191}
assert screen["support_sha256_sorted_xy_be24"]==prior["support_sha256_sorted_xy_be24"]
assert len(screen["overlap_counts"])==12
assert all(v==0 for v in screen["overlap_counts"].values())
assert screen["collision_relation_rank"]==0 and not screen["witnesses"]
assert sage["status"]=="LINE_SIGNED_SHORT_SAGE_INDEPENDENT_PASS"
assert sage["source_sha256"]==sha(HERE/"verify_line_signed_short_sage.py")
assert sage["input_sha256"]==sha(input_path)
assert sage["screen_sha256"]==sha(screen_path)
assert sage["runtime_sha256"]==sha(HERE/"type_i_line_signed_sage_runtime_info.json")
assert sage["support_counts"]==screen["support_counts"]
assert sage["overlap_counts"]==screen["overlap_counts"]
assert sage["support_sha256_sorted_xy_be24"]==screen["support_sha256_sorted_xy_be24"]
assert "LINE_SIGNED_SHORT_SAGE_INDEPENDENT_PASS" in (
    HERE/"type_i_line_signed_short_sage_verification.log").read_text()
assert "KeyError: 'support_sha256_sorted_xy_be24'" in (
    HERE/"type_i_line_signed_short_sage_verification.stale-producer.log").read_text()
print("LINE_SIGNED_SHORT_RECEIPT_PASS",json.dumps({
    "classes":141,"screen_cells":12,"overlaps":0,"rank":0},sort_keys=True))
