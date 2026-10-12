#!/usr/bin/env python3
"""Check CI metadata: routine PR jobs must have a proven budget <=25 minutes.

Manual-only, disabled and explicit label experiments keep their full budgets.
Unknown conditions are treated as applicable; unknown budgets fail closed.
This reads YAML/expressions only and never executes a workflow command.
"""
import ast
import hashlib
from pathlib import Path
import re
import sys

import yaml

UNKNOWN = object()
EVENTS = ("pull_request", "merge_group")
ROUTINE_PR_TYPES = {"opened", "synchronize", "reopened", "ready_for_review"}

# This one existing PR workflow is a historical proof artifact: its exact bytes
# are included in a checked Lean evidence seal. Keep the exception tied to the
# exact seal and workflow content, rather than allowing a workflow name or a
# proof-job label to exempt future work.
FROZEN_LEAN_WORKFLOW = ".github/workflows/weak-curve-lean.yml"
FROZEN_LEAN_JOB = "proofs"
FROZEN_LEAN_TIMEOUT = 60
FROZEN_LEAN_SEAL = "research/weak_curve_properties_lean_20261011/EVIDENCE_SHA256.tsv"
FROZEN_LEAN_SEAL_SHA256 = "118adb8dc0f96521b54b947989c4f3f84ecb25491de45563d71aef51efb6feea"
FROZEN_LEAN_PIN_PATH = "../../.github/workflows/weak-curve-lean.yml"
REPO_ROOT = Path(__file__).resolve().parents[2]


def frozen_lean_budget_exception(label, job_name, budget, event, root=REPO_ROOT):
    """Allow only the unchanged, evidence-sealed Lean proof job its 60m cap."""
    if (Path(label).as_posix() != FROZEN_LEAN_WORKFLOW
            or job_name != FROZEN_LEAN_JOB
            or budget != FROZEN_LEAN_TIMEOUT
            or event != "pull_request"):
        return False
    seal_path = root / FROZEN_LEAN_SEAL
    workflow_path = root / FROZEN_LEAN_WORKFLOW
    try:
        seal = seal_path.read_bytes()
        if hashlib.sha256(seal).hexdigest() != FROZEN_LEAN_SEAL_SHA256:
            return False
        pins = [line.split("\t", 1) for line in seal.decode("utf-8").splitlines()]
        pinned_hashes = [digest for digest, path in pins if path == FROZEN_LEAN_PIN_PATH]
        if len(pinned_hashes) != 1:
            return False
        return hashlib.sha256(workflow_path.read_bytes()).hexdigest() == pinned_hashes[0]
    except (OSError, UnicodeError, ValueError):
        return False


def expression(value, event, action=None):
    if not isinstance(value, str):
        return value
    source = value.strip()
    if source.startswith("${{") and source.endswith("}}"):
        source = source[3:-2].strip()
    # GitHub strings containing operators are outside this small evaluator.
    # Keep them unknown rather than guessing their meaning.
    for quoted in re.findall(r"'[^']*'|\"[^\"]*\"", source):
        if any(operator in quoted for operator in ("&&", "||", "!")):
            return UNKNOWN
    source = source.replace("&&", " and ").replace("||", " or ")
    source = re.sub(r"!(?!=)", " not ", source).strip()

    def visit(node):
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name):
            return {"true": True, "false": False, "null": None}.get(node.id, UNKNOWN)
        if isinstance(node, ast.Attribute):
            parts = []
            while isinstance(node, ast.Attribute):
                parts.insert(0, node.attr)
                node = node.value
            if isinstance(node, ast.Name):
                parts.insert(0, node.id)
            identity = ".".join(parts)
            if identity == "github.event_name":
                return event
            if identity == "github.event.action" and action is not None:
                return action
            return UNKNOWN
        if isinstance(node, ast.Call):
            return True if isinstance(node.func, ast.Name) and node.func.id == "always" and not node.args and not node.keywords else UNKNOWN
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
            operand = visit(node.operand)
            return UNKNOWN if operand is UNKNOWN else not operand
        if isinstance(node, ast.Compare) and len(node.ops) == 1:
            left, right = visit(node.left), visit(node.comparators[0])
            if left is UNKNOWN or right is UNKNOWN:
                return UNKNOWN
            if isinstance(left, str) and isinstance(right, str):
                left, right = left.casefold(), right.casefold()
            if isinstance(node.ops[0], ast.Eq):
                return left == right
            if isinstance(node.ops[0], ast.NotEq):
                return left != right
            return UNKNOWN
        if isinstance(node, ast.BoolOp):
            values = [visit(item) for item in node.values]
            unknown = False
            for item in values:
                if item is UNKNOWN:
                    unknown = True
                elif isinstance(node.op, ast.And) and not item:
                    return item
                elif isinstance(node.op, ast.Or) and item:
                    return UNKNOWN if unknown else item
            return UNKNOWN if unknown else values[-1]
        return UNKNOWN

    try:
        return visit(ast.parse("(" + source + ")", mode="eval").body)
    except (SyntaxError, ValueError):
        return UNKNOWN


def problems(document, label="workflow"):
    triggers = document.get("on", document.get(True, {}))
    if isinstance(triggers, str):
        triggers = {triggers: None}
    elif isinstance(triggers, list):
        triggers = dict.fromkeys(triggers)
    jobs = document.get("jobs", {})
    found = []
    for event in EVENTS:
        if event not in triggers:
            continue
        options = triggers[event] or {}
        types = options.get("types", []) if isinstance(options, dict) else []
        if event == "pull_request" and types and not ROUTINE_PR_TYPES.intersection(types):
            continue
        actions = sorted(ROUTINE_PR_TYPES.intersection(types)) if types else ["opened", "synchronize", "reopened"]
        if event == "merge_group":
            actions = ["checks_requested"]
        for action in actions:
            found.extend(job_problems(jobs, event, action, label))
    return found


def job_problems(jobs, event, action, label):
    applicable = {}

    def may_run(name):
        if name in applicable:
            return applicable[name]
        # Conservatively allow a malformed dependency cycle/missing job.
        applicable[name] = True
        job = jobs.get(name, {})
        condition = job.get("if", True)
        result = expression(condition, event, action)
        if result is not UNKNOWN and not result:
            applicable[name] = False
        else:
            needs = job.get("needs", [])
            needs = [needs] if isinstance(needs, str) else needs
            # Status functions can override default skipped-needs handling.
            overrides = isinstance(condition, str) and re.search(r"\b(always|failure|cancelled|success)\s*\(", condition)
            if not overrides and any(not may_run(parent) for parent in needs):
                applicable[name] = False
        return applicable[name]

    found = []
    for name, job in jobs.items():
        if not may_run(name):
            continue
        budget = expression(job.get("timeout-minutes", UNKNOWN), event, action)
        if isinstance(budget, bool) or not isinstance(budget, (int, float)) or not 0 < budget <= 25:
            if frozen_lean_budget_exception(label, name, budget, event):
                continue
            found.append(f"{label}: {name} can run on {event}/{action} without a proven 1–25 minute timeout")
    return found

def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".github/workflows")
    files = sorted(root.glob("*.yml")) + sorted(root.glob("*.yaml"))
    if not files:
        raise SystemExit(f"no workflows found in {root}")
    failures = []
    for path in files:
        failures.extend(problems(yaml.safe_load(path.read_text()), str(path)))
    for failure in failures:
        print(failure, file=sys.stderr)
    print(f"Checked routine PR/merge-group budgets in {len(files)} workflows: {len(failures)} violations")
    raise SystemExit(bool(failures))


if __name__ == "__main__":
    main()
