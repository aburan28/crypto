"""Configuration from the environment, shared by the CLI, worker and MCP server.

    TASKQ_REDIS_URL   redis://[:password@]host:6379/0   (required)
    TASKQ_NAMESPACE   key prefix, default `taskq`; one per deployment
    TASKQ_REPOS       path to a YAML/JSON map {alias: git url or local path}
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from .store import Store

DEFAULT_REPOS = {
    "crypto": "https://github.com/aburan28/crypto.git",
    "crypto-autoresearcher": "https://github.com/aburan28/crypto-autoresearcher.git",
    "cryptanalysis": "https://github.com/aburan28/cryptanalysis.git",
}


def store_from_env() -> Store:
    url = os.environ.get("TASKQ_REDIS_URL")
    if not url:
        raise SystemExit("TASKQ_REDIS_URL is not set")
    return Store.from_url(url, os.environ.get("TASKQ_NAMESPACE", "taskq"))


def load_repos(path: str | None = None) -> dict[str, str]:
    path = path or os.environ.get("TASKQ_REPOS")
    if not path:
        return dict(DEFAULT_REPOS)
    text = Path(path).read_text()
    if path.endswith((".yaml", ".yml")):
        import yaml
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    return {str(k): str(v) for k, v in (data or {}).items()}
