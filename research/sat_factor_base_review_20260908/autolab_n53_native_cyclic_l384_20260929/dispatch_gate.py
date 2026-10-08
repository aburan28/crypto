#!/usr/bin/env python3
"""Fail-closed, durable one-shot gate for the labeled hosted measurement."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import traceback
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
WORKFLOW = "n53-native-cyclic-l384.yml"
LABEL = "run-n53-native-cyclic-l384-once"


def api(path: str, token: str):
    request = Request(
        "https://api.github.com" + path,
        headers={"Accept": "application/vnd.github+json",
                 "Authorization": "Bearer " + token,
                 "X-GitHub-Api-Version": "2022-11-28"},
    )
    with urlopen(request, timeout=20) as response:
        return json.load(response)


def run(panel: Path):
    assert not panel.exists(), "refuse to reuse a panel path"
    panel.mkdir(parents=True)
    receipt_path = panel / "predispatch.json"
    receipt = {"schema": "n53_native_cyclic_l384_dispatch_gate_v1",
               "status": "CHECKING", "label": LABEL,
               "run_id": os.environ.get("GITHUB_RUN_ID"),
               "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT")}
    def save():
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    save()
    try:
        import check_protocol
        frozen = check_protocol.preflight(require_release=True)
        repo = os.environ["GITHUB_REPOSITORY"]
        run_id = int(os.environ["GITHUB_RUN_ID"])
        assert os.environ["GITHUB_RUN_ATTEMPT"] == "1", "rerun is not a new experiment"
        event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
        assert event["action"] == "labeled" and event["label"]["name"] == LABEL
        pr = event["pull_request"]
        number = pr["number"]
        head = pr["head"]["sha"]
        branch = pr["head"]["ref"]
        assert pr["head"]["repo"]["full_name"] == repo
        assert branch == frozen["expected_branch"]
        assert head == os.environ["KIC_NATIVE_L384_EXPECTED_HEAD"]
        live = api(f"/repos/{repo}/pulls/{number}", os.environ["GITHUB_TOKEN"])
        assert live["state"] == "open" and live["head"]["sha"] == head
        assert not (HERE / "evidence/archive_manifest.json").exists()
        prior_outcomes = []
        observed = 0
        for page in range(1, 11):
            runs = api(
                f"/repos/{repo}/actions/workflows/{WORKFLOW}/runs"
                f"?branch={branch}&event=pull_request&per_page=100&page={page}",
                os.environ["GITHUB_TOKEN"],
            )
            rows = runs["workflow_runs"]
            observed += len(rows)
            assert runs["total_count"] <= 1000, "one-shot history exceeds audit cap"
            for previous in rows:
                if previous["id"] == run_id or previous["head_branch"] != branch:
                    continue
                jobs = api(f"/repos/{repo}/actions/runs/{previous['id']}/jobs?per_page=100",
                           os.environ["GITHUB_TOKEN"])
                assert jobs["total_count"] <= 100
                for job in jobs["jobs"]:
                    if job["name"] == "outcome" and job["conclusion"] != "skipped":
                        prior_outcomes.append({"run_id": previous["id"],
                                               "job_id": job["id"],
                                               "status": job["status"],
                                               "conclusion": job["conclusion"]})
            if len(rows) < 100:
                break
        assert not prior_outcomes, f"prior outcome jobs exist: {prior_outcomes}"
        receipt.update({"status": "ADMITTED", "repo": repo, "pr_number": number,
                        "reviewed_head": head, "release_main_head": frozen["release_main_head"],
                        "prior_outcome_jobs": [], "workflow_runs_scanned": observed})
        save()
        print(json.dumps({"gate": "ADMITTED", "reviewed_head": head,
                          "prior_outcome_jobs": 0}, sort_keys=True))
    except Exception as exc:
        receipt.update({"status": "REFUSED",
                        "reason": f"{type(exc).__name__}: {exc}",
                        "traceback": traceback.format_exc()})
        save()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    args = parser.parse_args()
    run(args.panel)


if __name__ == "__main__":
    main()
