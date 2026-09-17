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

No type hints, camelCase identifiers (project convention).
"""
import subprocess
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
FLEET = (HERE / "fleet.sh").read_text()


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
        self.assertIn("FulfilledCapacity", roll)
        self.assertIn("${timeout}", roll)

    def test_roll_waits_against_a_per_batch_baseline_not_the_static_target(self):
        # Bugbot bc7697e0 (High): waiting for TargetCapacitySpecification
        # instead of the capacity this batch actually removed means an
        # already-short spot fleet -- one that never reaches that target --
        # times out on every batch and just keeps terminating instances,
        # draining the fleet instead of rolling it.
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertNotIn("TargetCapacitySpecification", roll)
        self.assertIn("before=$(capacityInt", roll)
        self.assertIn('"$fulfilled" -lt "$before"', roll)
        # The baseline must be read before terminate-instances runs.
        self.assertLess(roll.index("before=$(capacityInt"), roll.index("terminate-instances"))

    def test_roll_confirms_a_real_drop_before_trusting_a_recovered_reading(self):
        # Bugbot e7a26c6d (High): terminate-instances is async, so an
        # immediate FulfilledCapacity read can still show the pre-terminate
        # value. Treating that stale "already at $before" reading as
        # confirmation lets the next batch start before this one is
        # actually gone.
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        drop_wait = roll.split('if [ "$fulfilled" -lt "$before" ]; then break; fi', 1)
        self.assertEqual(len(drop_wait), 2, "no explicit wait for a drop below $before")
        # And only then a second wait for the climb back up.
        self.assertIn('while [ "$fulfilled" -lt "$before" ] && [ "$waited" -lt "$timeout" ]; do', roll)

    def test_roll_refuses_without_an_active_fleet(self):
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertIn('"$id" = None ]; then echo "no active fleet"', roll)


if __name__ == "__main__":
    unittest.main()
