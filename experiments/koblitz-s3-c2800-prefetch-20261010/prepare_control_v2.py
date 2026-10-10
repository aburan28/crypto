"""Freeze a direct-child RSS rerun after the macOS time wrapper failed."""

import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
OLD = HERE / "CONTROL_FREEZE.json"
NEW = HERE / "CONTROL_FREEZE_V2.json"
RUNNER = HERE / "run_control_v2.py"
FAILED = HERE / "local/first_wrapper_failure/receipt.json"


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    if NEW.exists():
        raise FileExistsError(NEW)
    original = OLD.read_bytes()
    failed = FAILED.read_bytes()
    row = json.loads(failed)
    assert row["freeze_sha256"] == sha(original)
    assert row["exit_code"] == 1 and row["raw_report_sha256"]
    assert row["unexpected_reference_differences"] == []
    assert row["peak_rss_bytes"] is None
    plan = json.loads(original)
    plan["kind"] = "s3_c2800_prefetch_t025_control_freeze_v2"
    plan["runner_sha256"] = sha(RUNNER.read_bytes())
    plan["parent_control_freeze_sha256"] = sha(original)
    plan["wrapper_failure_receipt_sha256"] = sha(failed)
    plan["resource_probe"] = "direct child plus resource.getrusage(RUSAGE_CHILDREN)"
    NEW.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"freeze_sha256": sha(NEW.read_bytes()),
                      "runner_sha256": plan["runner_sha256"]}, sort_keys=True))


if __name__ == "__main__":
    main()
