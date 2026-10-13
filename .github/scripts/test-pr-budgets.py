#!/usr/bin/env python3
"""Routing regressions for the CI metadata policy, using small YAML fixtures."""
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

spec = importlib.util.spec_from_file_location("budget_policy", Path(__file__).with_name("check-pr-budgets.py"))
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)


def workflow(job, event="pull_request", **jobs):
    return {"on": {event: None}, "jobs": {"check": job, **jobs}}


class BudgetPolicyTests(unittest.TestCase):
    def test_long_and_missing_automatic_budgets_fail(self):
        for budget in (None, 30, 75, 300):
            job = {} if budget is None else {"timeout-minutes": budget}
            self.assertTrue(policy.problems(workflow(job)))

    def test_short_budgets_pass(self):
        for budget in (1, 5, 12, 25):
            self.assertFalse(policy.problems(workflow({"timeout-minutes": budget})))

    def test_pr_full_event_expression_preserves_long_manual_lane(self):
        job = {"timeout-minutes": "${{ github.event_name == 'pull_request' && 25 || 300 }}"}
        self.assertFalse(policy.problems(workflow(job)))
        self.assertEqual(policy.expression(job["timeout-minutes"], "workflow_dispatch"), 300)

    def test_merge_group_requires_its_own_bounded_branch(self):
        job = {"timeout-minutes": "${{ github.event_name == 'pull_request' && 25 || 300 }}"}
        self.assertTrue(policy.problems(workflow(job, "merge_group")))
        job["timeout-minutes"] = "${{ (github.event_name == 'pull_request' || github.event_name == 'merge_group') && 25 || 300 }}"
        self.assertFalse(policy.problems(workflow(job, "merge_group")))

    def test_manual_and_disabled_jobs_keep_long_budgets(self):
        for condition in (False, "${{ false }}", "false && github.event_name == 'workflow_dispatch'",
                          "github.event_name == 'workflow_dispatch' && inputs.phase == 'full'",
                          "github.event_name != 'pull_request'"):
            self.assertFalse(policy.problems(workflow({"if": condition, "timeout-minutes": 300})))

    def test_unknown_or_conditions_cannot_hide_applicable_jobs(self):
        for condition in ("matrix.enabled", "github.event_name == 'workflow_dispatch' || matrix.enabled",
                          "github.event_name == 'workflow_dispatch' || github.event_name == 'pull_request'",
                          "unsupported.function()"):
            self.assertTrue(policy.problems(workflow({"if": condition, "timeout-minutes": 75})))

    def test_unknown_budgets_fail_closed(self):
        self.assertTrue(policy.problems(workflow({"timeout-minutes": "${{ matrix.timeout }}"})))

    def test_explicit_labeled_experiments_are_not_automatic_merge_checks(self):
        doc = workflow({"timeout-minutes": 300})
        doc["on"]["pull_request"] = {"types": ["labeled"]}
        self.assertFalse(policy.problems(doc))
        doc["on"]["pull_request"]["types"].append("synchronize")
        self.assertTrue(policy.problems(doc))

    def test_job_labels_do_not_exempt_automatic_or_branches(self):
        label_condition = "github.event.action == 'labeled' && github.event.label.name == 'run-once'"
        doc = workflow({"if": label_condition, "timeout-minutes": 300})
        self.assertFalse(policy.problems(doc))
        doc["jobs"]["check"]["if"] = "github.event.action == 'opened' || (" + label_condition + ")"
        self.assertTrue(policy.problems(doc))

    def test_multiline_conditions_preserve_manual_only_routing(self):
        doc = workflow({"if": "github.event_name == 'workflow_dispatch' &&\n inputs.phase == 'full'", "timeout-minutes": 300})
        self.assertFalse(policy.problems(doc))

    def test_github_string_comparisons_ignore_case(self):
        doc = workflow({"if": "github.event_name == 'PULL_REQUEST'", "timeout-minutes": 75})
        self.assertTrue(policy.problems(doc))

    def test_dependencies_inherit_manual_only_admission(self):
        doc = workflow({"if": "github.event_name == 'workflow_dispatch'", "timeout-minutes": 5},
                       producer={"needs": "check", "timeout-minutes": 300})
        self.assertFalse(policy.problems(doc))
        doc["jobs"]["producer"]["if"] = "always()"
        self.assertTrue(policy.problems(doc))

    def test_frozen_lean_exceptions_are_exact_and_seal_bound(self):
        for workflow_path, exception in policy.FROZEN_BUDGET_EXCEPTIONS.items():
            job = exception["job"]
            timeout = exception["timeout"]
            self.assertTrue(policy.frozen_proof_budget_exception(
                workflow_path, job, timeout, "pull_request"))
            self.assertFalse(policy.frozen_proof_budget_exception(
                workflow_path, "other", timeout, "pull_request"))
            self.assertFalse(policy.frozen_proof_budget_exception(
                workflow_path, job, timeout, "merge_group"))
            self.assertFalse(policy.frozen_proof_budget_exception(
                workflow_path, job, timeout + 1, "pull_request"))
            self.assertTrue(policy.problems(
                workflow({"timeout-minutes": 60}), workflow_path))

            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                seal = root / exception["seal"]
                workflow_file = root / workflow_path
                seal.parent.mkdir(parents=True)
                workflow_file.parent.mkdir(parents=True)
                shutil.copy2(policy.REPO_ROOT / exception["seal"], seal)
                shutil.copy2(policy.REPO_ROOT / workflow_path, workflow_file)
                self.assertTrue(policy.frozen_proof_budget_exception(
                    workflow_path, job, timeout, "pull_request", root=root))
                pinned_seal = seal.read_bytes()
                seal.write_bytes(pinned_seal + b"\n")
                self.assertFalse(policy.frozen_proof_budget_exception(
                    workflow_path, job, timeout, "pull_request", root=root))
                seal.write_bytes(pinned_seal)
                workflow_file.write_bytes(workflow_file.read_bytes() + b"\n# modified\n")
                self.assertFalse(policy.frozen_proof_budget_exception(
                    workflow_path, job, timeout, "pull_request", root=root))

    def test_metadata_is_not_executed(self):
        self.assertIs(policy.expression("__import__('os').system('false')", "pull_request"), policy.UNKNOWN)


if __name__ == "__main__":
    unittest.main()
