# Thin wrapper around the Rust rl_server binary.
#
#   GbrlEnv      — one subprocess, sequential rollouts.
#   VectorEnv    — N subprocesses, parallel rollouts. The trick is fire-and-forget
#                  writes followed by blocking reads: all N rl_servers compute
#                  their Buchberger step simultaneously while Python blocks on
#                  each readline in turn, so the total wall time per batch step
#                  is max(per-env step time), not sum. For cyclic-6+ where the
#                  Rust step is the bottleneck, this is ~N× faster.
#
# Optional: VectorEnv.step_partial supports a wall-clock timeout per step. Any
# proc that doesn't reply within `timeout_s` is killed and respawned; its slot
# in the result is None so the caller can finalize the trajectory.

import json
import select
import subprocess
from dataclasses import dataclass

import numpy as np


@dataclass
class StepResult:
    obs: dict
    reward: float
    done: bool


class GbrlEnv:
    def __init__(self, binary_path: str):
        self.proc = subprocess.Popen(
            [binary_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

    def _call(self, req: dict) -> dict:
        line = json.dumps(req) + "\n"
        self.proc.stdin.write(line)
        self.proc.stdin.flush()
        resp_line = self.proc.stdout.readline()
        if not resp_line:
            stderr = self.proc.stderr.read()
            raise RuntimeError(f"rl_server closed unexpectedly. stderr:\n{stderr}")
        resp = json.loads(resp_line)
        if not resp.get("ok"):
            raise RuntimeError(f"rl_server error: {resp.get('error')}")
        return resp

    def reset(self, instance: str) -> dict:
        return self._call({"method": "Reset", "instance": instance})["obs"]

    def step(self, pair_idx: int) -> StepResult:
        resp = self._call({"method": "Step", "pair_idx": pair_idx})
        return StepResult(obs=resp["obs"], reward=resp["reward"], done=resp["done"])

    def close(self):
        try:
            self.proc.stdin.close()
            self.proc.wait(timeout=2)
        except Exception:
            self.proc.kill()


class VectorEnv:
    """N parallel rl_server subprocesses. Step calls are fanned out write-all
    then read-all, so the N Rust workers compute in parallel.

    The active_mask argument to step_partial lets you advance only a subset of
    the envs in a batch — used when some envs have terminated mid-batch."""

    def __init__(self, binary_path: str, n: int):
        assert n >= 1
        self.n = n
        self.binary_path = binary_path
        self.procs = [self._spawn() for _ in range(n)]

    def _spawn(self) -> subprocess.Popen:
        return subprocess.Popen(
            [self.binary_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

    def _respawn(self, i: int):
        """Replace proc i with a fresh subprocess (after kill / hang)."""
        try:
            self.procs[i].kill()
            self.procs[i].wait(timeout=1)
        except Exception:
            pass
        self.procs[i] = self._spawn()

    def _send(self, i: int, req: dict):
        self.procs[i].stdin.write(json.dumps(req) + "\n")
        self.procs[i].stdin.flush()

    def _read(self, i: int) -> dict:
        line = self.procs[i].stdout.readline()
        if not line:
            err = self.procs[i].stderr.read()
            raise RuntimeError(f"VectorEnv: proc {i} closed unexpectedly. stderr:\n{err}")
        resp = json.loads(line)
        if not resp.get("ok"):
            raise RuntimeError(f"VectorEnv: proc {i} error: {resp.get('error')}")
        return resp

    def reset(self, instances) -> list:
        """Reset each subprocess to the corresponding instance.
        `instances` must have length self.n."""
        assert len(instances) == self.n
        for i, inst in enumerate(instances):
            self._send(i, {"method": "Reset", "instance": inst})
        return [self._read(i)["obs"] for i in range(self.n)]

    def reset_single(self, i: int, instance: str) -> dict:
        """Reset just proc i to a new instance. Used by async rollouts that
        immediately reuse a finished proc for the next trajectory."""
        self._send(i, {"method": "Reset", "instance": instance})
        return self._read(i)["obs"]

    def step_partial(self, actions, active_mask, timeout_s: float = None):
        """Send Step requests only for envs where active_mask[i] is True. Returns
        a list of length n where inactive slots are None, active slots are the
        full response dict from rl_server.

        If `timeout_s` is set, any proc that doesn't reply within that wall-clock
        budget is killed, respawned, and its slot is set to None — caller should
        treat that env as terminated and not assume the trajectory continues."""
        assert len(actions) == self.n and len(active_mask) == self.n
        # Phase 1: fire all writes (non-blocking buffered IO).
        for i in range(self.n):
            if active_mask[i]:
                self._send(i, {"method": "Step", "pair_idx": int(actions[i])})
        # Phase 2: collect responses. Each rl_server is computing in parallel
        # while we wait sequentially, so total wall = max single step time.
        out = [None] * self.n
        for i in range(self.n):
            if not active_mask[i]:
                continue
            if timeout_s is None:
                out[i] = self._read(i)
                continue
            # Poll stdout fd with select; respawn if no data within the budget.
            fd = self.procs[i].stdout
            ready, _, _ = select.select([fd], [], [], timeout_s)
            if not ready:
                self._respawn(i)
                out[i] = None
            else:
                try:
                    out[i] = self._read(i)
                except Exception:
                    self._respawn(i)
                    out[i] = None
        return out

    def close(self):
        for p in self.procs:
            try:
                p.stdin.close()
                p.wait(timeout=2)
            except Exception:
                p.kill()


# Per-feature scaling so raw integer degrees / term counts don't blow up the
# network. Picked from observed ranges on cyclic-n; tune per benchmark family.
PAIR_SCALE = np.array([
    10.0,  # lcm_degree
    10.0,  # sugar
    10.0,  # parent_i total_degree
    20.0,  # parent_i n_terms
    10.0,  # parent_i lm_degree
    10.0,  # parent_i sugar
    10.0,  # parent_j total_degree
    20.0,  # parent_j n_terms
    10.0,  # parent_j lm_degree
    10.0,  # parent_j sugar
], dtype=np.float32)

GLOBAL_SCALE = np.array([20.0, 50.0, 1.0], dtype=np.float32)


# Convert an Observation dict into (pair_feats[n_pairs, P], global_feats[G]).
# Returns None if the state has no critical pairs (terminal).
def featurize(obs: dict):
    pairs = obs["pairs"]
    polys = obs["polys"]
    if not pairs:
        return None

    pf_by_idx = {p["idx"]: p for p in polys}
    n_pairs = len(pairs)
    pair_feats = np.zeros((n_pairs, 10), dtype=np.float32)
    for k, pair in enumerate(pairs):
        pi = pf_by_idx[pair["i"]]
        pj = pf_by_idx[pair["j"]]
        pair_feats[k] = (
            pair["lcm_degree"],
            pair["sugar"],
            pi["total_degree"], pi["n_terms"], pi["lm_degree"], pi["sugar"],
            pj["total_degree"], pj["n_terms"], pj["lm_degree"], pj["sugar"],
        )
    pair_feats /= PAIR_SCALE

    global_feats = np.array([
        obs["basis_size"],
        n_pairs,
        np.log1p(obs["step_count"]),
    ], dtype=np.float32) / GLOBAL_SCALE

    return pair_feats, global_feats
