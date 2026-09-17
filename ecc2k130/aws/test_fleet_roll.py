"""Static contract for fleet.sh roll: the only way to deploy worker.py.

worker.py has no self-update path -- bootstrap.sh fetches it once, at
instance launch, and a running instance keeps executing whatever supervisor
code it booted with no matter how long `infra.sh sync` has since published a
fix (unlike the CUDA client, which workers poll campaign.json for and
self-restart onto). `roll` is what actually gets a worker.py fix onto an
already-running fleet: terminate instances in batches so `maintain` replaces
them and each re-bootstraps.  These checks pin that it terminates rather than
just describing, that it is genuinely batched (not everything at once, which
would drop the whole fleet's in-flight work at the same moment), that it
waits for each batch before starting the next, and that it is syntactically
sound bash.

The wait-loop tests below run the real `roll)` path against a mock `aws` so
they catch the two refill bugs static substring checks cannot: a stale
still-full `FulfilledCapacity` must not skip the wait, and the restore
target is the pre-batch fulfilled count, not `TotalTargetCapacity`.

No type hints, camelCase identifiers (project convention).
"""
import json
import os
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLEET = (HERE / "fleet.sh").read_text()

# Mock `aws` for roll: stateful EC2 Fleet view that can lag terminate-instances.
_MOCK_AWS = r"""#!/usr/bin/env python3
import json, os, sys
from pathlib import Path

statePath = Path(os.environ["ROLL_MOCK_STATE"])
logPath = Path(os.environ["ROLL_MOCK_LOG"])
state = json.loads(statePath.read_text())
args = sys.argv[1:]
logPath.write_text(logPath.read_text() + " ".join(args) + "\n")
joined = " ".join(args)

def save():
    statePath.write_text(json.dumps(state))

def idsText(ids):
    return "\t".join(ids)

def applyPending():
    pending = state.get("pending") or []
    if not pending:
        return
    gone = set(pending)
    state["instances"] = [i for i in state["instances"] if i not in gone]
    n = len(pending)
    # staleCapacity: terminate-instances has not yet moved FulfilledCapacity.
    if not state.get("staleCapacity"):
        state["fulfilled"] = max(0, int(state["fulfilled"]) - n)
    for extra in state.get("extras") or []:
        if extra not in state["instances"]:
            state["instances"].append(extra)
    if state.get("replace", True):
        nxt = int(state.get("nextId", 1))
        for _ in pending:
            state["instances"].append("i-new%03d" % nxt)
            nxt += 1
            if not state.get("staleCapacity"):
                state["fulfilled"] = int(state["fulfilled"]) + 1
        state["nextId"] = nxt
    state["pending"] = []

if len(args) >= 2 and args[0] == "ec2" and args[1] == "terminate-instances":
    ids = []
    if "--instance-ids" in args:
        i = args.index("--instance-ids") + 1
        while i < len(args) and not args[i].startswith("-"):
            ids.append(args[i])
            i += 1
    state.setdefault("terminates", []).append(ids)
    state["pending"] = ids
    state["staleLeft"] = int(state.get("stalePolls", 1))
    save()
    sys.exit(0)

if "ActiveInstances[].InstanceId" in joined:
    # Consume one stale wait-iteration after both describes have seen the
    # pre-terminate view (FulfilledCapacity is queried first in the loop).
    if state.get("mode") != "stale" and state.get("pending"):
        if int(state.get("staleLeft", 0)) > 0:
            state["staleLeft"] = int(state["staleLeft"]) - 1
            save()
        else:
            applyPending()
            save()
    sys.stdout.write(idsText(state["instances"]))
    sys.exit(0)

if "FulfilledCapacity" in joined:
    if state.get("mode") != "stale" and state.get("pending") and int(state.get("staleLeft", 0)) <= 0:
        applyPending()
        save()
    sys.stdout.write(str(state["fulfilled"]))
    sys.exit(0)

if "TotalTargetCapacity" in joined:
    sys.stdout.write(str(state.get("target", 8)))
    sys.exit(0)

sys.stderr.write("unexpected aws call: %s\n" % joined)
sys.exit(1)
"""

_MOCK_SLEEP = r"""#!/bin/sh
echo "$1" >> "$SLEEP_LOG"
"""


class FleetRollContract(unittest.TestCase):
    def test_syntax(self):
        subprocess.run(["bash", "-n", str(HERE / "fleet.sh")], check=True)

    def test_roll_is_documented(self):
        self.assertIn("./fleet.sh roll", FLEET)
        self.assertIn("ROLL_BATCH", FLEET)
        self.assertIn("ROLL_TIMEOUT_SECONDS", FLEET)

    def test_usage_fallback_covers_the_new_lines(self):
        # `*)` prints a fixed line range of this file as usage; it must not
        # go stale and clip roll's own doc line off the bottom.
        self.assertIn("sed -n '3,20p' \"$0\"", FLEET)
        head = FLEET.splitlines()[2:20]
        self.assertTrue(any("fleet.sh roll" in line for line in head))

    def test_roll_terminates_rather_than_just_reporting(self):
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertIn("aws ec2 terminate-instances", roll)

    def test_roll_batches_instead_of_terminating_everything_at_once(self):
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertIn("batch=${ROLL_BATCH:-4}", roll)
        # The chunking loop, not a single all-instance terminate call.
        self.assertIn('[ "$n" -lt "$batch" ]', roll)
        self.assertNotIn("--instance-ids $instances", roll)

    def test_roll_waits_for_the_fleet_to_refill_between_batches(self):
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        terminateAt = roll.index("aws ec2 terminate-instances")
        before = roll[:terminateAt]
        wait = roll[terminateAt:]
        # Restore target is snapshotted before terminate, from the live
        # fulfilled count — not TotalTargetCapacity, which a short Spot
        # fleet may never reach.
        self.assertIn("restore=", before)
        self.assertIn("FulfilledCapacity", before)
        self.assertIn("FulfilledCapacity", wait)
        self.assertIn('"$fulfilled" -ge "$restore"', wait)
        self.assertNotIn('"$fulfilled" -ge "$target"', wait)
        self.assertNotIn("TargetCapacitySpecification.TotalTargetCapacity", wait)
        # Stale FulfilledCapacity still equals the pre-terminate reading;
        # refill requires the terminated IDs to have left ActiveInstances.
        self.assertIn("ActiveInstances[].InstanceId", wait)
        # One unrelated new ID is not this batch's refill.
        self.assertIn('"$new" -ge "$n"', wait)
        self.assertNotIn('"$new" -eq 1', wait)
        # Nor is any 1-GPU capacity blip.
        self.assertIn("restore - n", wait)
        self.assertNotIn('"$fulfilled" -lt "$restore"', wait)
        self.assertIn("${timeout}", wait)

    def test_roll_refuses_without_an_active_fleet(self):
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertIn('"$id" = None ]; then echo "no active fleet"', roll)


class FleetRollWait(unittest.TestCase):
    def runRoll(self, state, batch=4, timeout=30):
        tmp = Path(tempfile.mkdtemp())
        try:
            script = tmp / "fleet.sh"
            script.write_text(FLEET)
            script.chmod(script.stat().st_mode | stat.S_IEXEC)
            aws = tmp / "aws"
            aws.write_text(_MOCK_AWS)
            aws.chmod(aws.stat().st_mode | stat.S_IEXEC)
            sleeper = tmp / "sleep"
            sleeper.write_text(_MOCK_SLEEP)
            sleeper.chmod(sleeper.stat().st_mode | stat.S_IEXEC)
            (tmp / ".fleet-id-us-west-2").write_text("fleet-test")
            statePath = tmp / "state.json"
            logPath = tmp / "aws.log"
            sleepLog = tmp / "sleep.log"
            statePath.write_text(json.dumps(state))
            logPath.write_text("")
            sleepLog.write_text("")
            env = os.environ.copy()
            env["PATH"] = str(tmp) + os.pathsep + env.get("PATH", "")
            env["ROLL_MOCK_STATE"] = str(statePath)
            env["ROLL_MOCK_LOG"] = str(logPath)
            env["SLEEP_LOG"] = str(sleepLog)
            env["ROLL_TIMEOUT_SECONDS"] = str(timeout)
            env["ROLL_BATCH"] = str(batch)
            env["AWS_DEFAULT_REGION"] = "us-west-2"
            proc = subprocess.run(
                ["bash", str(script), "roll"],
                cwd=str(tmp), env=env, capture_output=True, text=True)
            return proc, json.loads(statePath.read_text()), sleepLog.read_text(), logPath.read_text()
        finally:
            subprocess.run(["rm", "-rf", str(tmp)], check=False)

    def test_stale_capacity_does_not_skip_the_batch_wait(self):
        # Fleet accounting still shows a full fleet after terminate-instances.
        # Rolling the next batch immediately would SIGTERM everyone at once.
        proc, state, sleeps, _log = self.runRoll({
            "instances": ["i-1", "i-2", "i-3", "i-4", "i-5", "i-6"],
            "fulfilled": 6,
            "target": 6,
            "mode": "stale",
            "terminates": [],
        }, batch=4, timeout=30)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(state["terminates"]), 2)
        self.assertEqual(state["terminates"][0], ["i-1", "i-2", "i-3", "i-4"])
        self.assertEqual(state["terminates"][1], ["i-5", "i-6"])
        # timeout 30s / 15s poll: two sleeps per batch before the escape hatch.
        self.assertGreaterEqual(len(sleeps.split()), 4)

    def test_undercapacity_restores_prebatch_not_total_target(self):
        # 4 GPUs running against a requested 8. Replacements restore 4, never 8.
        # Waiting for TotalTargetCapacity would hit ROLL_TIMEOUT_SECONDS each
        # batch and drain the fleet; waiting for the pre-batch level does not.
        proc, state, sleeps, log = self.runRoll({
            "instances": ["i-1", "i-2", "i-3", "i-4"],
            "fulfilled": 4,
            "target": 8,
            "mode": "recover",
            "stalePolls": 1,
            "replace": True,
            "nextId": 1,
            "terminates": [],
        }, batch=2, timeout=1800)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(state["terminates"]), 2)
        self.assertNotIn("TotalTargetCapacity", log)
        self.assertLess(len(sleeps.split()), 20)
        self.assertGreaterEqual(len(sleeps.split()), 2)

    def test_unrelated_new_instance_does_not_end_the_batch_wait(self):
        # Terminated IDs have left ActiveInstances, FulfilledCapacity is
        # still the stale pre-batch reading, and one unrelated instance
        # appeared (shortfall fill, ReplaceUnhealthyInstances, or a
        # non-group interruption). That is not this batch's refill;
        # starting the next terminate now would overlap batches.
        proc, state, sleeps, _log = self.runRoll({
            "instances": ["i-1", "i-2", "i-3", "i-4", "i-5", "i-6"],
            "fulfilled": 6,
            "target": 6,
            "stalePolls": 1,
            "staleCapacity": True,
            "replace": False,
            "extras": ["i-extra"],
            "terminates": [],
        }, batch=4, timeout=30)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(state["terminates"]), 2)
        self.assertEqual(proc.stderr.count("moving on to the next batch"), 2)
        # timeout 30s / 15s poll: two sleeps per batch before the escape hatch.
        self.assertGreaterEqual(len(sleeps.split()), 4)


if __name__ == "__main__":
    unittest.main()
