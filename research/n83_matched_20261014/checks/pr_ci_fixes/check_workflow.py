#!/usr/bin/env python3
"""Check the CI configuration and its unchanged proof-source rejection gate."""
from pathlib import Path
import subprocess
import tempfile

import yaml


REPO = Path(__file__).resolve().parents[4]
WORKFLOW = REPO / ".github/workflows/n83-scalar-arithmetic-lean.yml"
AUDIT = REPO / "research/n83_speed_20261012/audit"
steps = yaml.safe_load(WORKFLOW.read_text())["jobs"]["proofs"]["steps"]
install_index = next(i for i, step in enumerate(steps)
                     if step.get("name") == "Install proof-source validation tooling")
guard_index = next(i for i, step in enumerate(steps)
                   if step.get("name") == "Reject incomplete or bypassed proof sources")
assert install_index < guard_index
assert "sudo apt-get update" in steps[install_index]["run"]
assert "sudo apt-get install -y ripgrep" in steps[install_index]["run"]
guard = steps[guard_index]["run"]


def check(directory, expected, env=None):
    result = subprocess.run(["/bin/bash", "-e", "-c", guard], cwd=directory,
                            env=env, text=True, capture_output=True)
    assert result.returncode == expected, result
    return result


check(AUDIT, 0)
print("The exact configured guard accepts the unchanged checked proof sources.")
with tempfile.TemporaryDirectory() as temporary:
    directory = Path(temporary)
    (directory / "AxiomAudit.lean").write_text("-- test fixture\n")
    for token in ("sorry", "admit", "sorryAx", "axiom", "native_decide", "#eval"):
        (directory / "ScalarArithmetic.lean").write_text(token + "\n")
        result = check(directory, 1)
        assert "Incomplete proof or unverified evaluation is not permitted." in result.stdout
        print(f"The exact configured guard rejects {token!r}.")
    (directory / "ScalarArithmetic.lean").write_text("-- test fixture\n")
    check(directory, 1, env={"PATH": str(directory)})
    print("A missing ripgrep executable fails closed; it cannot skip validation.")
print("Workflow YAML parses, installation precedes validation, and all gate checks pass.")
