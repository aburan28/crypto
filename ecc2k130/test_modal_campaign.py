#!/usr/bin/env python3
"""The campaign driver against a fake of the four SDK calls it makes."""
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import modal_campaign as mc  # noqa: E402


class FakeCalls:
    def __init__(self):
        self.spawned = []          # kwargs per spawn, in order
        self.results = {}          # call_id -> result dict, exception, or mc.Running
        self.cancelled = []
        self.n = 0

    def spawn(self, **kwargs):
        self.n += 1
        cid = "fc-%03d" % self.n
        self.spawned.append((cid, kwargs))
        self.results[cid] = mc.Running()
        return cid

    def poll(self, call_id):
        outcome = self.results[call_id]
        if isinstance(outcome, BaseException):
            raise outcome
        return outcome

    def cancel(self, call_id):
        self.cancelled.append(call_id)


class Clock:
    def __init__(self, t=1000.0):
        self.t = t

    def __call__(self):
        return self.t


def finished(dp=10, solved=None, cpu=None):
    return {"gpu": "RTX PRO 6000", "distinguishedPoints": dp, "rate": 14600.0,
            "iterations": 4 * 10 ** 12, "stopped": "deadline", "solved": solved, "cpu": cpu,
            "checkpoint": "/data/ckpt/curve131-run8000.ck"}


class DriverTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.state = os.path.join(self.tmp.name, "calls.json")
        self.calls = FakeCalls()
        self.clock = Clock()
        self.addCleanup(setattr, mc, "log", mc.log)
        self.logged = []
        mc.log = self.logged.append

    def driver(self, run_ids=(8000, 8001), hours=4.0, passes=0, **kwargs):
        base = dict(curve=131, packed=True, cpuThreads=32)
        base.update(kwargs)
        return mc.Driver(self.calls, run_ids, base, self.state, hours, passes,
                         grace_s=1800, now=self.clock)

    def test_first_tick_spawns_every_run_with_the_shape_and_its_own_id(self):
        d = self.driver()
        self.assertEqual(d.tick(), 2)
        ids = [kw["runId"] for _, kw in self.calls.spawned]
        self.assertEqual(ids, [8000, 8001])
        for _, kw in self.calls.spawned:
            self.assertEqual((kw["hours"], kw["packed"], kw["cpuThreads"]), (4.0, True, 32))
        saved = mc.load_state(self.state)
        self.assertEqual(saved["8000"]["call_id"], "fc-001")
        self.assertEqual(saved["8001"]["pass"], 1)

    def test_a_restart_reattaches_to_running_calls_instead_of_spawning_again(self):
        # The regression this driver exists for: a second container on the
        # same run id would fight the first over one checkpoint.
        self.driver().tick()
        again = self.driver()
        self.assertEqual(again.tick(), 2)
        self.assertEqual(len(self.calls.spawned), 2)

    def test_a_finished_pass_is_followed_by_the_next_at_once(self):
        d = self.driver()
        d.tick()
        self.calls.results["fc-001"] = finished(dp=1234)
        self.assertEqual(d.tick(), 2)
        self.assertEqual(len(self.calls.spawned), 3)
        cid, kw = self.calls.spawned[-1]
        self.assertEqual(kw["runId"], 8000)
        self.assertEqual(mc.load_state(self.state)["8000"], {"call_id": cid, "spawned_at": 1000.0, "pass": 2})
        self.assertTrue(any("1234 dp" in line for line in self.logged))

    def test_a_failed_or_unknown_call_is_held_without_a_second_container(self):
        d = self.driver()
        d.tick()
        self.calls.results["fc-002"] = RuntimeError("container lost")
        self.assertEqual(d.tick(), 1)
        self.assertEqual(len(self.calls.spawned), 2)
        self.assertTrue(mc.load_state(self.state)['8001']['attention'])
        self.assertEqual(d.run(poll_s=0), 1)
        self.assertTrue(any("failed" in line and "container lost" in line for line in self.logged))

    def test_a_call_past_its_deadline_is_cancelled_but_not_blindly_replaced(self):
        d = self.driver()
        d.tick()
        self.clock.t = 1000.0 + 4 * 3600 + 1799
        self.assertEqual(d.tick(), 2)
        self.assertEqual(self.calls.cancelled, [])
        self.clock.t = 1000.0 + 4 * 3600 + 1801
        self.assertEqual(d.tick(), 0)
        self.assertEqual(sorted(self.calls.cancelled), ["fc-001", "fc-002"])
        self.assertEqual(len(self.calls.spawned), 2)

    def test_cancellation_failure_never_spawns_an_overlapping_call(self):
        d = self.driver(run_ids=(8000,))
        d.tick()
        self.clock.t += 4 * 3600 + 1801
        def unavailable(call_id):
            raise RuntimeError('network error; old container may still be running')
        self.calls.cancel = unavailable
        self.assertEqual(d.tick(), 0)
        self.assertEqual(len(self.calls.spawned), 1)
        # Restarting the driver preserves the hold.
        self.assertEqual(self.driver(run_ids=(8000,)).run(poll_s=0), 1)
        self.assertEqual(len(self.calls.spawned), 1)

    def test_returned_error_or_missing_checkpoint_does_not_restart_seeds(self):
        for result in ({'error': 'build failed'}, {'stopped': 'deadline'},
                       {'returncode': 6, 'checkpoint': 'incompatible.ck'}, None):
            with self.subTest(result=result):
                if os.path.exists(self.state):
                    os.unlink(self.state)
                self.calls = FakeCalls()
                d = self.driver(run_ids=(8000,))
                d.tick()
                self.calls.results['fc-001'] = result
                self.assertEqual(d.tick(), 0)
                self.assertEqual(len(self.calls.spawned), 1)

    def test_a_solved_pass_stops_everything(self):
        d = self.driver()
        d.tick()
        self.calls.results["fc-001"] = finished(solved="k = 42")
        self.assertEqual(d.tick(), 0)
        self.assertEqual(d.solved["solved"], "k = 42")
        self.assertEqual(len(self.calls.spawned), 2)
        self.assertEqual(d.run(poll_s=0, sleep=lambda s: None), 0)

    def test_the_pass_budget_is_honoured(self):
        d = self.driver(run_ids=(8000,), passes=2)
        d.tick()
        self.calls.results["fc-001"] = finished()
        self.assertEqual(d.tick(), 1)                 # pass 2 spawned
        self.calls.results["fc-002"] = finished()
        self.assertEqual(d.tick(), 0)                 # no pass 3
        self.assertEqual(len(self.calls.spawned), 2)

    def test_summary_names_the_cpu_walker_and_a_skipped_one(self):
        line = mc.summarize(finished(cpu={"runId": 9000, "threads": 32, "build": "v4",
                                          "distinguishedPoints": 40, "rate": 400.0,
                                          "stopped": "deadline", "skipped": None}))
        self.assertIn("cpu run 9000: 32 threads v4, 40 dp, 400.0 M it/s", line)
        line = mc.summarize(finished(cpu={"runId": 9000, "threads": 32, "build": "v4",
                                          "distinguishedPoints": 0, "rate": 0,
                                          "stopped": "never started",
                                          "skipped": "host binary failed --test"}))
        self.assertIn("SKIPPED: host binary failed --test", line)
        self.assertIn("ERROR: build broke", mc.summarize({"error": "build broke"}))

    def test_cli_maps_the_shape_to_run_search_kwargs(self):
        args = mc.parse_args(["--run-id-base", "8000", "--count", "4", "--hours", "4",
                              "--packed", "--dp-weight", "32", "--cpu-threads", "32"])
        kw = mc.search_kwargs(args)
        self.assertEqual((kw["curve"], kw["packed"], kw["dpWeight"], kw["cpuThreads"],
                          kw["walksTarget"], kw["loadMax"]),
                         (131, True, 32, 32, 6160384, 2000000))
        # -1 leaves the image's ECC_CPU_THREADS in charge.
        kw = mc.search_kwargs(mc.parse_args(["--run-id-base", "8000"]))
        self.assertNotIn("cpuThreads", kw)


if __name__ == "__main__":
    raise SystemExit(unittest.main())
