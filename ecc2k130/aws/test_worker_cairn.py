"""The optional heartbeat to a cairn node's dashboard.

A fleet that walks ECC2K-130 for the campaign store can also be seen on a
cairn node's `/ui/task?id=...` page, which puts what the log paid a worker
beside what the worker says it is doing.  The second half is this heartbeat.
Three things are pinned here: the configuration comes from the environment
first and campaign.json second, and only when both the node and the
objective are named; the body is the shape the node's POST /progress
takes, built from counters this supervisor already has; and nothing the
node does -- down, slow, refusing, missing the route -- reaches the
campaign loop, which is the one property that matters, since the lease
heartbeat right before it is what keeps the slot.

No network: urlopen is replaced.  No type hints, camelCase identifiers
(project convention).
"""
import io
import json
import os
import shutil
import tempfile
import unittest
import urllib.error
from unittest import mock

import worker
from worker import (CAIRN_CLIENT, CAIRN_DEFAULT_TRAIL_BITS, cairnConfig,
                    cairnHeartbeat, postCairnHeartbeat)


class Config(unittest.TestCase):
    def test_needs_both_a_node_and_an_objective(self):
        self.assertIsNone(cairnConfig({}, env={}))
        self.assertIsNone(cairnConfig({"cairnNode": "http://n:8080"}, env={}))
        self.assertIsNone(cairnConfig({}, env={"ECC_CAIRN_OBJECTIVE": "sha256:o"}))
        conf = cairnConfig({"cairnNode": "http://n:8080/", "cairnObjective": "sha256:o"}, env={})
        self.assertEqual(conf["node"], "http://n:8080", "no trailing slash: the path is appended")
        self.assertEqual(conf["objective"], "sha256:o")
        self.assertEqual(conf["worker"], "")
        self.assertEqual(conf["trailBits"], CAIRN_DEFAULT_TRAIL_BITS)

    def test_the_environment_wins_over_campaign_json(self):
        cfg = {"cairnNode": "http://fleet:8080", "cairnObjective": "sha256:fleet",
               "cairnWorker": "fleet", "cairnTrailBits": 12}
        env = {"ECC_CAIRN_NODE": "http://mine:9090", "ECC_CAIRN_WORKER": "me"}
        conf = cairnConfig(cfg, env=env)
        self.assertEqual(conf["node"], "http://mine:9090")
        self.assertEqual(conf["objective"], "sha256:fleet", "unset in the environment, so the file's")
        self.assertEqual(conf["worker"], "me")
        self.assertEqual(conf["trailBits"], 12)

    def test_a_bad_trail_bits_falls_back_rather_than_raising(self):
        conf = cairnConfig({"cairnNode": "http://n", "cairnObjective": "o", "cairnTrailBits": "x"}, env={})
        self.assertEqual(conf["trailBits"], CAIRN_DEFAULT_TRAIL_BITS)


class Body(unittest.TestCase):
    conf = {"node": "http://n:8080", "objective": "sha256:o", "worker": "", "trailBits": 16}

    def test_a_run_that_has_printed_reports_its_own_counters(self):
        last = {"rate": 14.003e9, "iters": 123456789012, "dp": 4321, "stored": 4300, "dropped": 2}
        body = cairnHeartbeat(self.conf, 7, {"dpUploaded": 4000, "walks": 385024}, last,
                              spoolBytes=72 * 300, recordBytes=72, device="NVIDIA RTX PRO 6000", lanes=385024)
        self.assertEqual(body, {
            "objective_id": "sha256:o",
            "worker": "slot-00007",
            "steps": 123456789012,
            "trails": 4321,
            "units_submitted": 4000,
            "units_pending": 300,
            "trail_bits": 16,
            "client": CAIRN_CLIENT,
            "steps_per_second": 14003000000,
            "device": "NVIDIA RTX PRO 6000",
            "lanes": 385024,
        })

    def test_a_fresh_run_says_it_is_here_and_nothing_more(self):
        body = cairnHeartbeat(dict(self.conf, worker="alice"), 7, {}, None, 0, 72, None, 0)
        self.assertEqual(body["worker"], "alice", "the pseudonym the orbits are submitted under")
        self.assertEqual(body["steps"], 0)
        self.assertEqual(body["trails"], 0)
        self.assertEqual(body["units_submitted"], 0)
        self.assertEqual(body["units_pending"], 0)
        for absent in ("steps_per_second", "device", "lanes"):
            self.assertNotIn(absent, body, "unknown is left out, not zeroed")

    def test_the_body_is_json_the_node_can_take(self):
        body = cairnHeartbeat(self.conf, 1, {}, {"rate": 1e6, "iters": 5, "dp": 1, "stored": 1, "dropped": 0},
                              0, 32, "cpu", 1)
        text = json.dumps(body)
        self.assertNotIn(".", json.dumps(body["steps_per_second"]), "integers only: cairn has no float")
        self.assertEqual(json.loads(text), body)


class FakeResponse(io.BytesIO):
    def __init__(self, status, payload):
        super().__init__(payload)
        self.status = status

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class Posting(unittest.TestCase):
    conf = {"node": "http://n:8080", "objective": "sha256:o", "worker": "w", "trailBits": 16}

    def test_posts_json_to_the_progress_route(self):
        seen = {}

        def fakeUrlopen(request, timeout=None):
            seen["url"] = request.full_url
            seen["method"] = request.get_method()
            seen["type"] = request.get_header("Content-type")
            seen["body"] = json.loads(request.data.decode())
            seen["timeout"] = timeout
            return FakeResponse(202, b'{"recorded":true,"status":"live"}')

        with mock.patch.object(worker.urllib.request, "urlopen", fakeUrlopen):
            status, text = postCairnHeartbeat(self.conf, {"objective_id": "sha256:o", "worker": "w", "steps": 1})
        self.assertEqual(status, 202)
        self.assertIn("recorded", text)
        self.assertEqual(seen["url"], "http://n:8080/progress")
        self.assertEqual(seen["method"], "POST")
        self.assertEqual(seen["type"], "application/json")
        self.assertEqual(seen["body"]["worker"], "w")
        self.assertEqual(seen["timeout"], worker.CAIRN_HEARTBEAT_TIMEOUT)

    def test_a_refusal_is_returned_not_raised(self):
        def refuse(request, timeout=None):
            raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {},
                                         io.BytesIO(b'{"error":"no such objective in this log"}'))

        with mock.patch.object(worker.urllib.request, "urlopen", refuse):
            status, text = postCairnHeartbeat(self.conf, {})
        self.assertEqual(status, 404)
        self.assertIn("no such objective", text)


class Supervisor(unittest.TestCase):
    """`Worker.cairnBeat` on a supervisor shaped by hand: no S3, no GPU."""

    def supervisor(self, cfg):
        w = worker.Worker.__new__(worker.Worker)
        w.cfg = cfg
        w.state = {"dpUploaded": 12, "walks": 8}
        w.gpuName = "test-gpu"
        # `dpPath` hangs off the work directory; an empty one means no dp
        # file yet, which `dpStride` answers with the v1 stride.
        w.work = tempfile.mkdtemp(prefix="cairn-beat-")
        self.addCleanup(shutil.rmtree, w.work, True)
        w.spoolBytes = lambda: 0
        w.cairnFailures = 0
        w.cairnOff = False
        return w

    def test_does_nothing_without_a_configuration(self):
        w = self.supervisor({})
        with mock.patch.object(worker, "postCairnHeartbeat") as post, \
                mock.patch.dict(os.environ, {}, clear=True):
            w.cairnBeat(3, None)
        post.assert_not_called()

    def test_posts_and_never_raises_whatever_the_node_does(self):
        cfg = {"cairnNode": "http://n:8080", "cairnObjective": "sha256:o"}
        w = self.supervisor(cfg)
        # A refusal first, then nine dead connections, then the node is back:
        # eleven ticks, ten of them failures, one success.
        outcomes = [(503, "busy")] + [Exception("connection refused")] * 9 + [(202, '{"recorded":true}')]
        calls = []

        def post(conf, body, timeout=None):
            calls.append(body)
            outcome = outcomes.pop(0)
            if isinstance(outcome, Exception):
                raise outcome
            return outcome

        logged = []
        with mock.patch.object(worker, "postCairnHeartbeat", post), \
                mock.patch.object(worker, "log", logged.append), \
                mock.patch.dict(os.environ, {}, clear=True):
            for _ in range(11):
                w.cairnBeat(3, {"rate": 2e9, "iters": 10, "dp": 1, "stored": 1, "dropped": 0})
        self.assertEqual(len(calls), 11, "every tick posts; a failure never stops the next")
        self.assertEqual(calls[0]["worker"], "slot-00003")
        self.assertEqual(calls[0]["units_submitted"], 12)
        self.assertEqual(calls[0]["lanes"], 8)
        # The first failure and every tenth are logged, not all ten: the
        # refusal is the first, the tenth dead connection the second.
        self.assertEqual(len(logged), 2, logged)
        self.assertIn("refused (503, 1 so far)", logged[0])
        self.assertIn("failed (10 so far", logged[1])
        self.assertEqual(w.cairnFailures, 0, "a 202 clears the count")
        self.assertFalse(w.cairnOff)

    def test_a_node_without_the_route_is_told_once_and_left_alone(self):
        cfg = {"cairnNode": "http://old:8080", "cairnObjective": "sha256:o"}
        w = self.supervisor(cfg)
        calls = []

        def post(conf, body, timeout=None):
            calls.append(body)
            return 404, '{"error":"no such path"}'

        logged = []
        with mock.patch.object(worker, "postCairnHeartbeat", post), \
                mock.patch.object(worker, "log", logged.append), \
                mock.patch.dict(os.environ, {}, clear=True):
            for _ in range(5):
                w.cairnBeat(3, None)
        self.assertEqual(len(calls), 1)
        self.assertTrue(w.cairnOff)
        self.assertTrue(any("no POST /progress" in line for line in logged), logged)

    def test_an_unknown_objective_is_a_refusal_that_keeps_trying(self):
        # The node is new enough; the objective is not in its log yet. That
        # can change (the coordinator posts it), so this is retried.
        cfg = {"cairnNode": "http://n:8080", "cairnObjective": "sha256:o"}
        w = self.supervisor(cfg)
        with mock.patch.object(worker, "postCairnHeartbeat",
                               lambda conf, body, timeout=None: (404, '{"error":"no such objective in this log"}')), \
                mock.patch.object(worker, "log", lambda line: None), \
                mock.patch.dict(os.environ, {}, clear=True):
            w.cairnBeat(3, None)
            w.cairnBeat(3, None)
        self.assertFalse(w.cairnOff)
        self.assertEqual(w.cairnFailures, 2)


if __name__ == "__main__":
    unittest.main()
