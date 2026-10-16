#!/usr/bin/env python3
"""Confirm an immutable GitHub selector artifact before fresh-point work."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import urllib.request

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPOSITORY = "aburan28/crypto"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_receipt(artifact, artifact_id, action_digest, run_id, attempt,
                  commit, selection, selector_host):
    normalized = action_digest.removeprefix("sha256:").lower()
    expected_name = f"n83-tier1-selector-{run_id}-{attempt}"
    workflow = artifact.get("workflow_run") or {}
    if (artifact.get("id") != artifact_id
            or artifact.get("name") != expected_name
            or artifact.get("expired") is not False
            or not isinstance(artifact.get("size_in_bytes"), int)
            or artifact["size_in_bytes"] <= 0
            or len(normalized) != 64
            or any(character not in "0123456789abcdef"
                   for character in normalized)
            or artifact.get("digest") != "sha256:" + normalized
            or workflow.get("id") != run_id
            or workflow.get("head_sha") != commit
            or selection.get("screen_commit") != commit
            or selector_host.get("git_sha") != commit
            or len(selector_host.get("cpu_affinity", [])) != 2
            or selector_host.get("cpu_envelope") != "github_two_core"):
        raise ValueError("uploaded GitHub selector artifact differs from frozen run")
    return {
        "schema": "n83.tier1-github-selector-seal/v1",
        "repository": REPOSITORY,
        "run_id": run_id, "run_attempt": attempt, "head_sha": commit,
        "artifact_id": artifact_id, "artifact_name": expected_name,
        "artifact_digest": "sha256:" + normalized,
        "artifact_size_bytes": artifact["size_in_bytes"],
        "artifact_url": artifact.get("archive_download_url"),
        "selection_sha256": digest(HERE / "selection.json"),
        "candidate_manifest_sha256": selection["candidate_manifest_sha256"],
        "selector_host_sha256": selection["screen_host_sha256"],
        "selector_analysis_sha256": selection["screen_analysis_sha256"],
        "selector_cpus": selector_host["cpu_affinity"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-id", required=True, type=int)
    parser.add_argument("--artifact-digest", required=True)
    parser.add_argument("--selector-results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    token = os.environ.get("GH_TOKEN")
    run = os.environ.get("GITHUB_RUN_ID", "")
    attempt = os.environ.get("GITHUB_RUN_ATTEMPT", "")
    if (not token or not run.isdecimal() or not attempt.isdecimal()
            or os.environ.get("GITHUB_REPOSITORY") != REPOSITORY):
        raise ValueError("GitHub token or run identity unavailable")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                     cwd=ROOT, text=True).strip()
    if os.environ.get("GITHUB_SHA") != commit:
        raise ValueError("GitHub run source commit differs")
    result = args.selector_results.resolve()
    selection = json.loads((HERE / "selection.json").read_text())
    selector_host = json.loads((result / "host.json").read_text())
    candidate_name = Path(selection["candidate_manifest"]).name
    sealed_copy = result.parent
    if (digest(sealed_copy / "selection.json") != digest(
                HERE / "selection.json")
            or digest(sealed_copy / candidate_name) != selection[
                "candidate_manifest_sha256"]):
        raise ValueError("uploaded selector or candidate copy differs")
    if (selection["screen_analysis_sha256"] != digest(result / "analysis.json")
            or selection["screen_source_sha256"] != digest(result / "source.json")
            or selection["screen_host_sha256"] != digest(result / "host.json")):
        raise ValueError("selector analysis, source or host changed before seal")
    request = urllib.request.Request(
        f"https://api.github.com/repos/{REPOSITORY}/actions/artifacts/"
        f"{args.artifact_id}",
        headers={"Accept": "application/vnd.github+json",
                 "Authorization": "Bearer " + token,
                 "X-GitHub-Api-Version": "2022-11-28"})
    with urllib.request.urlopen(request, timeout=30) as response:
        artifact = json.load(response)
    receipt = build_receipt(artifact, args.artifact_id,
                            args.artifact_digest, int(run), int(attempt),
                            commit, selection, selector_host)
    blob = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    for path in (HERE / "github_artifact_seal.json", args.output.resolve()):
        if path.exists() and path.read_text() != blob:
            raise ValueError(f"GitHub selector seal changed: {path}")
        path.write_text(blob)
    print(json.dumps({key: receipt[key] for key in
                     ("artifact_id", "artifact_digest", "selection_sha256")},
                     sort_keys=True))


if __name__ == "__main__":
    main()
