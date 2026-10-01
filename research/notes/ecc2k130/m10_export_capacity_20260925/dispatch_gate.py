#!/usr/bin/env python3
"""One-shot reviewed-head admission for the hosted n131 representation-size run."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import traceback
from urllib.request import Request, urlopen

import run

HERE = Path(__file__).resolve().parent
WORKFLOW = "ecc2k130-m10-capacity-once.yml"
LABEL = "run-ecc2k130-m10-capacity-once"
BRANCH = "codex/n131-m10-capacity-attempt-20260929"
REVIEWER_ASSOCIATIONS = {"OWNER", "MEMBER", "COLLABORATOR"}


def api(path: str, token: str):
    request = Request("https://api.github.com" + path, headers={
        "Accept": "application/vnd.github+json",
        "Authorization": "Bearer " + token,
        "X-GitHub-Api-Version": "2022-11-28",
    })
    with urlopen(request, timeout=20) as response:
        return json.load(response)


def latest_reviews(rows: list[dict]) -> dict[str, dict]:
    latest = {}
    for row in sorted(rows, key=lambda item: item.get("submitted_at") or ""):
        if row["state"] in {"APPROVED", "CHANGES_REQUESTED", "DISMISSED"}:
            latest[row["user"]["login"]] = row
    return latest


def check_prior_capacity(repo: str, branch: str, current_run: int, token: str) -> int:
    observed = 0
    for page in range(1, 11):
        runs = api(
            f"/repos/{repo}/actions/workflows/{WORKFLOW}/runs"
            f"?branch={branch}&event=pull_request&per_page=100&page={page}", token)
        rows = runs["workflow_runs"]
        observed += len(rows)
        assert runs["total_count"] <= 1000, "one-shot history exceeds audit cap"
        for previous in rows:
            if previous["id"] == current_run or previous["head_branch"] != branch:
                continue
            jobs = api(f"/repos/{repo}/actions/runs/{previous['id']}/jobs?per_page=100", token)
            assert jobs["total_count"] <= 100
            for job in jobs["jobs"]:
                if job["name"] == "capacity" and job["conclusion"] != "skipped":
                    raise RuntimeError(f"prior capacity job exists: {previous['id']}/{job['id']}")
        if len(rows) < 100:
            break
    return observed


def admit(out: Path) -> dict:
    assert not out.exists(), "refusing to overwrite a gate receipt directory"
    out.mkdir(parents=True)
    receipt_path = out / "predispatch.json"
    receipt = {"schema": "ecc2k130-m10-capacity-dispatch-v1",
               "utc": datetime.now(timezone.utc).isoformat(),
               "decision": "REFUSED", "measurement_children": 0,
               "run_id": os.environ.get("GITHUB_RUN_ID"),
               "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT")}
    def save() -> None:
        receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    save()
    try:
        token = os.environ["GITHUB_TOKEN"]
        repo = os.environ["GITHUB_REPOSITORY"]
        current_run = int(os.environ["GITHUB_RUN_ID"])
        assert os.environ["GITHUB_RUN_ATTEMPT"] == "1", "workflow rerun is not a new attempt"
        event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
        assert event["action"] == "labeled" and event["label"]["name"] == LABEL
        pr = event["pull_request"]
        pr_number = int(pr["number"])
        assert pr_number == int(event["number"])
        assert pr["head"]["ref"] == BRANCH
        assert pr["head"]["repo"]["full_name"] == repo
        head = pr["head"]["sha"]
        assert head == subprocess.check_output(["git", "rev-parse", "HEAD"],
                                               cwd=run.ROOT, text=True).strip()
        live = api(f"/repos/{repo}/pulls/{pr_number}", token)
        assert live["state"] == "open" and not live["draft"]
        assert live["number"] == pr_number
        assert live["head"]["sha"] == head and live["head"]["ref"] == BRANCH
        reviews = api(f"/repos/{repo}/pulls/{pr_number}/reviews?per_page=100", token)
        assert len(reviews) < 100, "review list needs pagination"
        latest = latest_reviews(reviews)
        assert not any(row["state"] == "CHANGES_REQUESTED" for row in latest.values())
        approved = [row for row in latest.values()
                    if row["state"] == "APPROVED" and row["commit_id"] == head
                    and row["user"]["id"] != live["user"]["id"]
                    and row["author_association"] in REVIEWER_ASSOCIATIONS]
        assert approved, "independent exact-head approval is missing"
        observed = check_prior_capacity(repo, BRANCH, current_run, token)
        frozen = json.loads((HERE / "FROZEN.json").read_text())
        release = run.release_gate(frozen)
        assert release["checkout_head"] == head
        receipt.update({"decision": "ADMITTED", "repo": repo,
                        "pr_number": pr_number, "reviewed_head": head,
                        "reviewer_logins": sorted(row["user"]["login"] for row in approved),
                        "prior_capacity_jobs": 0,
                        "workflow_runs_scanned": observed,
                        "release": release})
        save()
        return receipt
    except BaseException as exc:
        receipt["error_type"] = type(exc).__name__
        receipt["error"] = str(exc)
        receipt["traceback"] = traceback.format_exc()
        save()
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = admit(args.out.resolve())
    print(json.dumps({"decision": result["decision"],
                      "reviewed_head": result["reviewed_head"]}, sort_keys=True))


if __name__ == "__main__":
    main()
