#!/usr/bin/env python3
"""Unauthenticated + env-based access probe for Modal and RunPod (no paid launch)."""
from __future__ import annotations

import json
import os
import subprocess
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parent / "access-probe.json"


def head(url: str, timeout: float = 10.0) -> dict:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return {"ok": True, "status": resp.status, "url": url}
    except urllib.error.HTTPError as e:
        return {"ok": False, "status": e.code, "url": url, "reason": str(e.reason)}
    except Exception as e:  # noqa: BLE001 — diagnostic probe
        return {"ok": False, "status": None, "url": url, "error_type": type(e).__name__, "error": str(e)}


def main() -> None:
    env_modal = {
        "MODAL_TOKEN_ID_set": bool(os.environ.get("MODAL_TOKEN_ID")),
        "MODAL_TOKEN_SECRET_set": bool(os.environ.get("MODAL_TOKEN_SECRET")),
    }
    env_runpod = {"RUNPOD_API_KEY_set": bool(os.environ.get("RUNPOD_API_KEY"))}
    pip = {}
    for pkg in ("modal", "runpod"):
        p = subprocess.run(
            ["python3", "-c", f"import {pkg}; print(getattr({pkg}, '__version__', 'ok'))"],
            capture_output=True,
            text=True,
        )
        pip[pkg] = {
            "import_ok": p.returncode == 0,
            "stdout": p.stdout.strip(),
            "stderr": p.stderr.strip()[:500],
        }
    payload = {
        "schema": "isogeny_rho_gpu_access_probe/v1",
        "utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "gpu_execution_confirmed": False,
        "env": {"modal": env_modal, "runpod": env_runpod},
        "pip": pip,
        "http": {
            "modal_api": head("https://api.modal.com"),
            "runpod_graphql": head("https://api.runpod.io/graphql"),
        },
        "prior_repo_receipt": {
            "path": "research/notes/ecc2k130/tau_adic_arithmetic_20260929/gpu/results/runpod-access-01.json",
            "status": "http_error",
            "http_status": 403,
            "note": "prior authenticated RunPod inventory failed Forbidden",
        },
        "decision": (
            "No Modal/RunPod credentials in this environment; do not launch a paid GPU job. "
            "Leaf-curve GPU rho is unimplemented in ecc2k130/; automorphism floor already "
            "prices leaf rho √131 worse than E0 before any kernel constant."
        ),
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
