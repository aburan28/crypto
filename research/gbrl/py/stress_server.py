# Stress-test rl_server: random walk on a long-running instance with periodic
# progress prints. Designed to surface the cyclic-5 hang from training.

import json
import os
import random
import subprocess
import sys
import time


def main():
    binary = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "target", "release", "rl_server"))
    instance = sys.argv[1] if len(sys.argv) > 1 else "cyclic-5"
    max_steps = int(sys.argv[2]) if len(sys.argv) > 2 else 50_000

    proc = subprocess.Popen(
        [binary],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )

    def call(req):
        t0 = time.time()
        proc.stdin.write(json.dumps(req) + "\n")
        proc.stdin.flush()
        line = proc.stdout.readline()
        dt = time.time() - t0
        if not line:
            err = proc.stderr.read()
            raise RuntimeError(f"server closed. stderr:\n{err}")
        return json.loads(line), dt

    rng = random.Random(0)
    t_total = time.time()
    print(f"reset {instance}", flush=True)
    resp, dt = call({"method": "Reset", "instance": instance})
    obs = resp["obs"]
    print(f"  reset: dt={dt*1000:.0f}ms basis={obs['basis_size']} pairs={len(obs['pairs'])}", flush=True)

    max_dt = 0.0
    max_obs_bytes = 0
    for step in range(max_steps):
        pairs = obs["pairs"]
        if not pairs:
            print(f"  terminated at step {step} (no pairs)", flush=True)
            break
        idx = rng.randrange(len(pairs))
        resp, dt = call({"method": "Step", "pair_idx": idx})
        obs = resp["obs"]
        max_dt = max(max_dt, dt)
        obs_bytes = len(json.dumps(obs))
        max_obs_bytes = max(max_obs_bytes, obs_bytes)
        if step % 500 == 0 or dt > 0.5:
            print(f"  step={step:5d} dt={dt*1000:6.1f}ms basis={obs['basis_size']:4d} "
                  f"pairs={len(pairs):4d} obs_bytes={obs_bytes}", flush=True)
        if obs["done"]:
            print(f"  done at step {step+1}", flush=True)
            break

    print(f"total wall: {time.time()-t_total:.1f}s  max step dt: {max_dt*1000:.0f}ms  "
          f"max obs bytes: {max_obs_bytes}", flush=True)
    proc.stdin.close()
    proc.wait(timeout=5)


if __name__ == "__main__":
    main()
