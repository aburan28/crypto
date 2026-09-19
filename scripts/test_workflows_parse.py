#!/usr/bin/env python3
"""Every workflow file must parse as YAML.

A workflow GitHub cannot parse does not fail loudly in the way a failing step
does: it produces a run with no jobs, named after the file path rather than the
workflow, and the checks inside it simply stop happening. That is how
`.github/workflows/polynomial-reuse.yml` spent four days red on every push to
main while its redis-cache correctness tests -- the whole point of the file --
ran not once. The cause was one character:

    - run: cargo test --release --features redis-cache --lib koblitz_groebner::

A plain YAML scalar may not end in a colon, so the trailing `::` of a module
filter turns the line into a mapping and the file into a parse error. It is
invisible in review, it does not resemble a broken pipeline, and nothing else
in 85 workflow files would report it.

This test is deliberately narrower than a workflow linter. It asserts only
that each file parses and has the two keys a workflow cannot work without,
because that is the failure mode that hides.
"""

import pathlib
import unittest

try:
    import yaml
except ImportError:  # pragma: no cover - exercised only where PyYAML is absent
    yaml = None

WORKFLOWS = pathlib.Path(__file__).resolve().parent.parent / ".github" / "workflows"


@unittest.skipIf(yaml is None, "PyYAML not installed")
class WorkflowsParse(unittest.TestCase):
    def files(self):
        found = sorted(WORKFLOWS.glob("*.yml")) + sorted(WORKFLOWS.glob("*.yaml"))
        self.assertTrue(found, "no workflow files found at %s" % WORKFLOWS)
        return found

    def test_every_workflow_parses_as_yaml(self):
        broken = []
        for path in self.files():
            try:
                yaml.safe_load(path.read_text(encoding="utf-8"))
            except yaml.YAMLError as err:
                mark = getattr(err, "problem_mark", None)
                broken.append("%s: %s%s" % (
                    path.name,
                    getattr(err, "problem", err),
                    " (line %d)" % (mark.line + 1) if mark else ""))
        self.assertEqual(broken, [], "unparseable workflow files: " + "; ".join(broken))

    def test_every_workflow_declares_a_trigger_and_a_job(self):
        for path in self.files():
            doc = yaml.safe_load(path.read_text(encoding="utf-8"))
            with self.subTest(workflow=path.name):
                self.assertIsInstance(doc, dict, "workflow is not a mapping")
                # PyYAML resolves the bare key `on` to the boolean True.
                self.assertTrue(doc.get("on") or doc.get(True), "no `on:` triggers")
                self.assertTrue(doc.get("jobs"), "no `jobs:`")

    def test_a_trailing_colon_in_a_run_step_is_caught(self):
        """The check discriminates: this is the exact text that broke main."""
        with self.assertRaises(yaml.YAMLError):
            yaml.safe_load(
                "jobs:\n  j:\n    steps:\n"
                "      - run: cargo test --lib koblitz_groebner::\n")
        yaml.safe_load(
            "jobs:\n  j:\n    steps:\n"
            '      - run: "cargo test --lib koblitz_groebner::"\n')


if __name__ == "__main__":
    unittest.main()
