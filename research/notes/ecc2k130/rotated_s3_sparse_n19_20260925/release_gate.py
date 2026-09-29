#!/usr/bin/env python3
"""Refuse a repeat or stale labeled first measurement of the frozen n19 gate."""
from __future__ import annotations

import json
import os
import subprocess
import time
import urllib.request
from pathlib import Path

import ci_replay
import run


def api(path: str):
    token = os.environ["GITHUB_TOKEN"]
    repository = os.environ["GITHUB_REPOSITORY"]
    request = urllib.request.Request(
        f"https://api.github.com/repos/{repository}/{path}",
        headers={"Accept": "application/vnd.github+json",
                 "Authorization": f"Bearer {token}",
                 "X-GitHub-Api-Version": "2022-11-28"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def main() -> None:
    frozen = ci_replay.check_freeze()
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    number = frozen["release_pr_number"]
    label = frozen["release_label"]
    assert os.environ["GITHUB_EVENT_NAME"] == "pull_request"
    assert os.environ["GITHUB_RUN_ATTEMPT"] == "1"
    assert os.environ["GITHUB_REPOSITORY"] == "aburan28/crypto"
    assert event["action"] == "labeled" and event["label"]["name"] == label
    assert event["number"] == number
    event_pr = event["pull_request"]
    assert event_pr["head"]["repo"]["full_name"] == "aburan28/crypto"
    assert event_pr["base"]["ref"] == "main"
    event_head = event_pr["head"]["sha"]
    assert subprocess.check_output(["git", "rev-parse", "HEAD"],
                                   text=True).strip() == event_head
    live_pr = api(f"pulls/{number}")
    assert live_pr["head"]["sha"] == event_head
    assert live_pr["base"]["ref"] == "main"
    assert live_pr["draft"] and live_pr["state"] == "open"
    assert label in {row["name"] for row in live_pr["labels"]}
    labeled = []
    for _ in range(6):
        labeled = []
        page = 1
        while True:
            rows = api(f"issues/{number}/timeline?per_page=100&page={page}")
            labeled.extend(row for row in rows if row.get("event") == "labeled"
                           and row.get("label", {}).get("name") == label)
            if len(rows) < 100:
                break
            page += 1
        if labeled:
            break
        time.sleep(5)
    assert len(labeled) == 1, "measurement label must have exactly one lifetime application"
    run.verify_parent_merged(frozen)
    run.require_linux_proc()
    print(json.dumps({"decision": "FIRST_RELEASE_APPROVED", "pr": number,
                      "head_sha": event_head, "label_event_id": labeled[0]["id"],
                      "run_id": int(os.environ["GITHUB_RUN_ID"]),
                      "run_attempt": int(os.environ["GITHUB_RUN_ATTEMPT"]),
                      "repository": os.environ["GITHUB_REPOSITORY"],
                      "freeze_sha256": ci_replay.sha(ci_replay.HERE / "FROZEN.json")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
