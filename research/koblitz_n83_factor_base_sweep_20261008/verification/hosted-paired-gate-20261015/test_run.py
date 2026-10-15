"""No-launch controls for the hosted pair's immutable inputs and cgroup gate."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from unittest import TestCase, main, mock


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("n83_hosted_pair_run", HERE / "run.py")
assert spec is not None and spec.loader is not None
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class HostedPairControls(TestCase):
    def test_wrong_frozen_plan_is_rejected_before_archive_replay(self) -> None:
        with mock.patch.object(runner, "sha256", return_value="0" * 64), \
             mock.patch.object(runner.preflight, "make_plan") as replay:
            with self.assertRaisesRegex(ValueError, "frozen paired plan bytes changed"):
                runner.checked_plan(Path("/unused"))
            replay.assert_not_called()

    def test_wrong_cgroup_is_rejected_before_solver_launch(self) -> None:
        result = mock.Mock(stdout="4294967296\n1073741824\n", stderr="")
        with mock.patch.object(runner.subprocess, "run", return_value=result) as call:
            with self.assertRaisesRegex(ValueError, "4 GiB and zero swap"):
                runner.cgroup_probe("sha256:" + "a" * 64)
        command = call.call_args.args[0]
        for fragment in ("--network", "none", "--cpus", "1", "--memory", "4096m",
                         "--memory-swap", "4096m"):
            self.assertIn(fragment, command)

    def test_non_linux_host_is_rejected_before_archive_io(self) -> None:
        with mock.patch.object(runner.platform, "system", return_value="Darwin"), \
             mock.patch.object(runner, "checked_plan") as plan:
            with self.assertRaisesRegex(ValueError, "requires Linux arm64"):
                runner.run(HERE, Path("/private/tmp/new-n83-hosted-output"),
                           Path("/private/tmp/new-n83-hosted-scratch"))
            plan.assert_not_called()


if __name__ == "__main__":
    main()
