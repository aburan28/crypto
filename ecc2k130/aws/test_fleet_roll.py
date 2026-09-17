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
        self.assertIn("activeInstanceIds", roll)
        self.assertIn("${timeout}", roll)

    def test_roll_waits_against_a_per_batch_baseline_not_the_static_target(self):
        # Bugbot bc7697e0 (High): waiting for TargetCapacitySpecification
        # instead of the count this batch actually removed means an
        # already-short spot fleet -- one that never reaches that target --
        # times out on every batch and just keeps terminating instances,
        # draining the fleet instead of rolling it.
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertNotIn("TargetCapacitySpecification", roll)
        self.assertIn("capacityBefore=$(activeGpuCapacity", roll)
        self.assertIn('"$capacity" -ge "$capacityBefore"', roll)
        # The baseline must be read before terminate-instances runs.
        self.assertLess(roll.index("capacityBefore=$(activeGpuCapacity"), roll.index("terminate-instances"))

    def test_roll_confirms_this_batchs_own_ids_are_gone_not_an_aggregate_number(self):
        # Bugbot e7a26c6d (High, first fix) and 2db76048 (High, on the fix
        # itself): FulfilledCapacity lags terminate-instances and unrelated
        # fleet churn can move it, so a still-full or "already recovered"
        # aggregate reading doesn't prove THIS batch is gone. Membership
        # (describe-fleet-instances) is the only thing that does.
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        # A comment may still explain why FulfilledCapacity is avoided; the
        # query itself must not be.
        self.assertNotIn("Fleets[0].FulfilledCapacity", roll)
        # Checks each terminated id's presence in the current active set,
        # not just a count, before declaring the batch gone.
        self.assertIn("for member in $group", roll)
        self.assertIn('case "$active" in *" $member "*) gone=0 ;; esac', roll)
        # Only after confirming departure does it wait for GPU-weighted
        # membership (not FulfilledCapacity, not instance count) to climb
        # back to the per-batch baseline.
        drop_then_recover = roll.split('if [ "$gone" -eq 1 ]; then break; fi', 1)
        self.assertEqual(len(drop_then_recover), 2)
        self.assertIn('if [ "$gone" -eq 1 ]; then', drop_then_recover[1])
        self.assertIn('"$capacity" -ge "$capacityBefore"', drop_then_recover[1])
        self.assertNotIn("wc -w", drop_then_recover[1])

    def test_roll_recovery_weights_by_gpu_not_instance_count(self):
        # Bugbot 36e4c416 (High): waiting for the active-instance count to
        # recover treats refill as instance-count, but the fleet's maintain
        # target and WeightedCapacity are GPUs. After terminating large
        # instances, a handful of smaller replacements restore the count
        # while most of the lost GPU capacity is still missing, and the
        # next batch then drains the campaign.
        start = FLEET.index("activeGpuCapacity()")
        helper = FLEET[start:FLEET.index("\ncapacityInt")]
        self.assertIn("gpusOf", helper)
        self.assertIn("ActiveInstances[].InstanceType", helper)
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertNotIn("wc -w", roll)
        self.assertIn("capacityBefore=$(activeGpuCapacity", roll)
        self.assertIn('"$capacity" -ge "$capacityBefore"', roll)

    def test_roll_propagates_a_failed_capacity_query_instead_of_reading_it_as_zero(self):
        # Bugbot 9ab64fed (High): a bare `types=$(aws ...)` inside a
        # function is not enough -- bash does not inherit `set -e` into a
        # command substitution by default, so activeGpuCapacity's own
        # internal errexit does nothing when it is itself invoked as
        # $(activeGpuCapacity ...) one level up. A failing describe-fleet-
        # instances call would silently produce total=0, and 0 always
        # satisfies "capacity has recovered", so the next batch would
        # proceed without the previous one's replacements ever landing.
        # The fix must check the inner command explicitly instead of
        # relying on inherited errexit.
        start = FLEET.index("activeGpuCapacity()")
        helper = FLEET[start:FLEET.index("\ncapacityInt")]
        self.assertIn("types=$(aws", helper)
        self.assertIn("|| return 1", helper)

    def test_roll_normalizes_the_tab_separated_id_list_before_matching(self):
        # aws --output text separates list entries with tabs, not spaces;
        # a bare `*" $member "*` substring check against that raw text only
        # ever matches the first and last id (the ones next to the literal
        # spaces this script adds), silently never matching any id in the
        # middle of the list.
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertIn("tr '\\t' ' '", roll)

    def test_roll_refuses_without_an_active_fleet(self):
        roll = FLEET.split("roll)", 1)[1].split("\n*)", 1)[0]
        self.assertIn('"$id" = None ]; then echo "no active fleet"', roll)


if __name__ == "__main__":
    unittest.main()
