#!/usr/bin/env python3
"""Keep the campaign's Modal runs walking, whether or not anyone is watching.

`modal run modal_app.py::fanout` owns an *ephemeral* app: when the local client
disconnects the app stops, and with --detach Modal keeps only the last spawned
function alive. On 2026-09-21 four campaign runners were terminated at 15:23Z
that way while the client kept drawing its spinner for four hours.

A *deployed* app has no such client. Its functions are spawned as calls that
run to completion on their own, and a call can be found again by id from any
process. This driver:

  * spawns one runSearch per run id on the deployed app and writes the call
    ids to a state file;
  * on every start re-attaches to the calls in that file that are still
    running, instead of spawning a second container onto the same checkpoint;
  * when a call returns (the pass reached its deadline), spawns the next pass
    for that run id at once, so the checkpoint is resumed within a minute;
  * holds a failed, unobservable or overdue call for recovery. An uncertain
    cancellation is not proof that the old container stopped, so it must
    never authorize a second container on those seeds.

    modal deploy modal_app.py                 # once per code change, same env as run.sh
    python3 modal_campaign.py --run-id-base 8000 --count 4 --hours 4 --packed ...

The shape flags are the ones ::search takes, and they must be the same every
pass: the checkpoint header records curve, run id, worker count and batch,
and the client refuses a checkpoint whose shape differs.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

APP_NAME = "ecc2k130"
FUNCTION = "runSearch"
DEFAULT_STATE = os.path.expanduser("~/.ecc2k130-modal-campaign/calls.json")
POLL_S = 60
# A pass that has not returned this long after its deadline is a lost
# container, not a slow checkpoint: the client stops itself at `hours` and
# gets ten minutes to checkpoint.
GRACE_S = 1800


def log(msg):
    sys.stderr.write(time.strftime("%Y-%m-%dT%H:%M:%SZ ", time.gmtime()) + msg + "\n")
    sys.stderr.flush()


class Running(Exception):
    """The call has not returned yet."""


class ModalCalls:
    """The SDK surface the driver uses, thin enough to fake in tests."""

    def __init__(self, app_name=APP_NAME, function=FUNCTION, environment=None):
        import modal
        self.modal = modal
        self.fn = modal.Function.from_name(app_name, function, environment_name=environment)

    def spawn(self, **kwargs):
        return self.fn.spawn(**kwargs).object_id

    def poll(self, call_id):
        """The call's result, or raise Running; any other exception is the
        call's own failure."""
        import modal
        call = modal.FunctionCall.from_id(call_id)
        try:
            return call.get(timeout=0)
        except modal.exception.TimeoutError:
            raise Running()
        except TimeoutError:
            raise Running()

    def cancel(self, call_id):
        import modal
        modal.FunctionCall.from_id(call_id).cancel(terminate_containers=True)


def load_state(path):
    if not os.path.isfile(path):
        return {}
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def save_state(path, state):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(state, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def summarize(result):
    """One line of what a pass did, from runSearch's return dict."""
    if not isinstance(result, dict):
        return repr(result)[:200]
    if result.get("error"):
        return "ERROR: " + str(result["error"])[-300:]
    line = ("%s: %s dp, %.1f B it/s, %s iters, stopped=%s"
            % (result.get("gpu", "?"), result.get("distinguishedPoints"),
               float(result.get("rate") or 0) / 1e3, result.get("iterations"),
               result.get("stopped")))
    cpu = result.get("cpu")
    if cpu:
        line += ("; cpu run %s: %s threads %s, %s dp, %.1f M it/s, stopped=%s"
                 % (cpu.get("runId"), cpu.get("threads"), cpu.get("build"),
                    cpu.get("distinguishedPoints"), float(cpu.get("rate") or 0),
                    cpu.get("stopped")))
        if cpu.get("skipped"):
            line += " SKIPPED: " + str(cpu["skipped"])[:200]
    if result.get("solved"):
        line += "  SOLVED: " + str(result["solved"])
    return line


class Driver:
    def __init__(self, calls, run_ids, kwargs, state_path, hours, passes=0,
                 grace_s=GRACE_S, now=time.time):
        self.calls = calls
        self.run_ids = [int(r) for r in run_ids]
        self.kwargs = dict(kwargs)
        self.state_path = state_path
        self.hours = float(hours)
        self.passes = int(passes)
        self.grace_s = grace_s
        self.now = now
        self.state = load_state(state_path)
        self.solved = None
        self.attention = None

    def hold(self, run_id, reason):
        entry = self.entry(run_id)
        entry["attention"] = reason
        save_state(self.state_path, self.state)
        self.attention = reason
        log("run %d held for operator recovery: %s" % (run_id, reason))

    def entry(self, run_id):
        return self.state.get(str(run_id))

    def spawn(self, run_id, reason):
        entry = self.entry(run_id) or {}
        pass_no = int(entry.get("pass") or 0) + 1
        if self.passes and pass_no > self.passes:
            log("run %d: %d passes done; not spawning another" % (run_id, self.passes))
            return False
        call_id = self.calls.spawn(runId=run_id, hours=self.hours, **self.kwargs)
        self.state[str(run_id)] = {"call_id": call_id, "spawned_at": self.now(),
                                   "pass": pass_no}
        save_state(self.state_path, self.state)
        log("run %d: pass %d spawned as %s (%s)" % (run_id, pass_no, call_id, reason))
        return True

    def tick(self):
        """Look at every run once. Returns the number of runs still walking."""
        walking = 0
        for run_id in self.run_ids:
            entry = self.entry(run_id)
            if entry and entry.get("attention"):
                self.attention = entry["attention"]
                continue
            if not entry or not entry.get("call_id"):
                if self.spawn(run_id, "no call on record"):
                    walking += 1
                continue
            try:
                result = self.calls.poll(entry["call_id"])
            except Running:
                age = self.now() - float(entry.get("spawned_at") or 0)
                if age > self.hours * 3600 + self.grace_s:
                    log("run %d: call %s is %.0f s past its deadline; cancelling"
                        % (run_id, entry["call_id"], age - self.hours * 3600))
                    try:
                        self.calls.cancel(entry["call_id"])
                    except Exception as exc:
                        log("run %d: cancel failed: %s: %s" % (run_id, type(exc).__name__, exc))
                    self.hold(run_id, "deadline exceeded; confirm the old container stopped and its checkpoint before replacement")
                else:
                    walking += 1
                continue
            except Exception as exc:
                log("run %d: pass %s failed: %s: %s" % (run_id, entry.get("pass"),
                                                       type(exc).__name__, str(exc)[:300]))
                self.hold(run_id, "call result failed or is unknown; inspect the existing call before replacement")
                continue
            log("run %d: pass %s finished -- %s" % (run_id, entry.get("pass"), summarize(result)))
            if isinstance(result, dict) and result.get("solved"):
                self.solved = result
                return 0
            if (not isinstance(result, dict) or result.get("error")
                    or result.get("returncode", 0) != 0 or not result.get("checkpoint")):
                self.hold(run_id, "pass returned without a successful checkpoint; refusing to restart seeds")
                continue
            if self.spawn(run_id, "next pass"):
                walking += 1
        return walking

    def run(self, poll_s=POLL_S, sleep=time.sleep):
        while True:
            walking = self.tick()
            if self.attention:
                log("operator recovery required; no replacement calls were spawned for held runs")
                return 1
            if self.solved:
                log("SOLVED: %s" % self.solved.get("solved"))
                return 0
            if walking == 0:
                log("nothing left to walk")
                return 0
            sleep(poll_s)


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-id-base", type=int, required=True)
    ap.add_argument("--count", type=int, default=1)
    ap.add_argument("--hours", type=float, default=4.0)
    ap.add_argument("--passes", type=int, default=0, help="0 = until solved or stopped")
    ap.add_argument("--curve", type=int, default=131)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--threads", type=int, default=256)
    ap.add_argument("--leaf", type=int, default=0)
    ap.add_argument("--walks", type=int, default=6160384)
    ap.add_argument("--dp-weight", type=int, default=-1)
    ap.add_argument("--load-max", type=int, default=2000000)
    ap.add_argument("--verify", type=int, default=0)
    ap.add_argument("--checkpoint-every", type=int, default=60)
    ap.add_argument("--packed", action="store_true")
    ap.add_argument("--off-campaign", action="store_true")
    ap.add_argument("--cpu-threads", type=int, default=-1,
                    help="-1: the deployed image's ECC_CPU_THREADS; 0: GPU only")
    ap.add_argument("--app", default=APP_NAME)
    ap.add_argument("--environment", default=None)
    ap.add_argument("--state", default=os.environ.get("ECC_MODAL_CAMPAIGN_STATE", DEFAULT_STATE))
    ap.add_argument("--poll", type=float, default=POLL_S)
    return ap.parse_args(argv)


def search_kwargs(args):
    kwargs = dict(curve=args.curve, batch=args.batch, threads=args.threads, leaf=args.leaf,
                  dpWeight=args.dp_weight, walksTarget=args.walks, loadMax=args.load_max,
                  packed=args.packed, verify=args.verify,
                  checkpointEvery=args.checkpoint_every, offCampaign=args.off_campaign)
    if args.cpu_threads >= 0:
        kwargs["cpuThreads"] = args.cpu_threads
    return kwargs


def main(argv=None):
    args = parse_args(argv)
    # The same rules modal_app.py applies at launch, checked here before any
    # container is rented: modal_app imports `modal`, which this host has.
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import modal_app
    run_ids = list(range(args.run_id_base, args.run_id_base + args.count))
    with_cpu = args.cpu_threads > 0 or (args.cpu_threads < 0 and modal_app.CPU_THREADS > 0)
    if modal_app.isCampaignRun(args.curve, args.off_campaign):
        for run_id in run_ids:
            modal_app.checkCampaignRunId(run_id, withCpu=with_cpu)
        if args.dp_weight >= 0:
            modal_app.campaignDpWeightFor(args.dp_weight)
    calls = ModalCalls(args.app, FUNCTION, args.environment)
    driver = Driver(calls, run_ids, search_kwargs(args), args.state, args.hours, args.passes)
    log("driving runs %s on deployed app %s: %.2f h passes, state %s"
        % (run_ids, args.app, args.hours, args.state))
    return driver.run(poll_s=args.poll)


if __name__ == "__main__":
    sys.exit(main())
