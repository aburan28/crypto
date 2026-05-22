# Minimal PPO trainer over the Buchberger env. Sequential rollouts (one env),
# one trajectory per update. Designed as a learning-signal smoke test, not a
# production trainer — but the architecture (pointer policy, variable action
# space, value baseline, GAE) is the right starting point.
#
#     python3 train.py --instances cyclic-4 --episodes 200

import argparse
import math
import os
import random
from collections import deque

import numpy as np
import torch
import torch.nn.functional as F
from torch.distributions import Categorical

from env import GbrlEnv, VectorEnv, featurize
from model import PointerPolicy


def rollout_async_parallel(
    vec_env: VectorEnv,
    policy: PointerPolicy,
    instance_pool,
    max_steps: int,
    num_trajs: int,
    device,
    step_timeout_s=None,
    reward_mode: str = "step",
):
    """Async rollouts: keep all N procs busy until `num_trajs` trajectories are
    collected. When one proc's rollout ends (terminal, max_steps, or step
    timeout), it's immediately reset to a fresh random instance from `pool` and
    starts the next rollout — no waiting for slow stragglers.

    Returns (entries, metrics) where:
      entries: list of (instance, traj, steps, status) of length num_trajs;
               status ∈ {"completed", "truncated", "timeout"}.
      metrics: dict with counts of each status."""
    import random as _rng
    n = vec_env.n
    instances = [_rng.choice(instance_pool) for _ in range(n)]
    obses = vec_env.reset(instances)
    trajs = [[] for _ in range(n)]
    insts = list(instances)
    steps_count = [0] * n
    # For reward_mode="arith_ops": track per-env cumulative arith_ops between
    # steps so reward[k] = -(arith_ops_after_step_k - arith_ops_after_step_{k-1}).
    prev_ops = [int(obs.get("arith_ops", 0)) for obs in obses]
    completed = []  # list of (instance, traj, steps, status)

    def finalize_and_reset(i, status):
        completed.append((insts[i], trajs[i], steps_count[i], status))
        if len(completed) >= num_trajs:
            return
        new_inst = _rng.choice(instance_pool)
        insts[i] = new_inst
        steps_count[i] = 0
        trajs[i] = []
        obses[i] = vec_env.reset_single(i, new_inst)
        prev_ops[i] = int(obses[i].get("arith_ops", 0))

    while len(completed) < num_trajs:
        actions = [0] * n
        step_mask = [False] * n
        for i in range(n):
            if len(completed) >= num_trajs:
                break
            # If current state is terminal or rollout cap hit, finalize and reset.
            f = featurize(obses[i])
            if f is None or steps_count[i] >= max_steps:
                status = "completed" if f is None else "truncated"
                finalize_and_reset(i, status)
                if len(completed) >= num_trajs:
                    break
                f = featurize(obses[i])
                if f is None:
                    # Trivial instance (no pairs at reset). Skip this proc tick.
                    continue
            pair_feats, global_feats = f
            pf = torch.from_numpy(pair_feats).to(device)
            gf = torch.from_numpy(global_feats).to(device)
            with torch.no_grad():
                logits, value = policy(pf, gf)
                dist = Categorical(logits=logits)
                action = dist.sample()
                logp = dist.log_prob(action)
            a = int(action.item())
            trajs[i].append({
                "pair_feats": pair_feats,
                "global_feats": global_feats,
                "action": a,
                "logp": float(logp.item()),
                "value": float(value.item()),
                "reward": 0.0,
            })
            actions[i] = a
            step_mask[i] = True

        if not any(step_mask):
            break

        results = vec_env.step_partial(actions, step_mask, timeout_s=step_timeout_s)
        for i in range(n):
            if not step_mask[i]:
                continue
            r = results[i]
            if r is None:
                # Step timeout: respawn already happened in vec_env. Pop the
                # un-rewarded sample, then finalize whatever we have.
                if trajs[i]:
                    trajs[i].pop()
                finalize_and_reset(i, "timeout")
            else:
                if reward_mode == "arith_ops":
                    new_ops = int(r["obs"].get("arith_ops", 0))
                    delta = new_ops - prev_ops[i]
                    # Penalize work done this step (more ops → more negative reward).
                    trajs[i][-1]["reward"] = float(-delta)
                    prev_ops[i] = new_ops
                else:
                    trajs[i][-1]["reward"] = r["reward"]
                obses[i] = r["obs"]
                steps_count[i] += 1

    completed = completed[:num_trajs]
    metrics = {
        "n_completed": sum(1 for _, _, _, s in completed if s == "completed"),
        "n_truncated": sum(1 for _, _, _, s in completed if s == "truncated"),
        "n_timeout":   sum(1 for _, _, _, s in completed if s == "timeout"),
        "completed_lens": [steps for _, _, steps, s in completed if s == "completed"],
        "all_lens": [steps for _, _, steps, _ in completed],
    }
    return completed, metrics


def gae(traj, gamma: float, lam: float):
    n = len(traj)
    values = [t["value"] for t in traj] + [0.0]
    rewards = [t["reward"] for t in traj]
    advs = [0.0] * n
    g = 0.0
    for t in reversed(range(n)):
        delta = rewards[t] + gamma * values[t + 1] - values[t]
        g = delta + gamma * lam * g
        advs[t] = g
    returns = [advs[t] + values[t] for t in range(n)]
    return advs, returns


def ppo_update_batched(policy, optimizer, batch, device,
                        epochs: int, clip: float, vc: float, ec: float):
    """Batched PPO update: pad all batch elements to max(n_pairs), run a single
    forward per epoch with masking. Replaces the per-step Python loop in the
    update phase — large speedup when batches are big."""
    B = len(batch)
    P_max = max(b["pair_feats"].shape[0] for b in batch)
    feat_dim = batch[0]["pair_feats"].shape[1]
    global_dim = batch[0]["global_feats"].shape[0]

    pair_feats_pad = np.zeros((B, P_max, feat_dim), dtype=np.float32)
    pair_mask = np.zeros((B, P_max), dtype=bool)
    global_feats_stack = np.zeros((B, global_dim), dtype=np.float32)
    for i, b in enumerate(batch):
        p = b["pair_feats"]
        pair_feats_pad[i, :p.shape[0]] = p
        pair_mask[i, :p.shape[0]] = True
        global_feats_stack[i] = b["global_feats"]

    pf_t = torch.from_numpy(pair_feats_pad).to(device)
    gf_t = torch.from_numpy(global_feats_stack).to(device)
    mask_t = torch.from_numpy(pair_mask).to(device)

    advs_t = torch.tensor([b["adv"] for b in batch], dtype=torch.float32, device=device)
    returns_t = torch.tensor([b["return"] for b in batch], dtype=torch.float32, device=device)
    old_logps = torch.tensor([b["logp"] for b in batch], dtype=torch.float32, device=device)
    actions = torch.tensor([b["action"] for b in batch], dtype=torch.long, device=device)

    last_pl = last_vl = last_ent = last_kl = last_logit_max = 0.0
    for _ in range(epochs):
        logits, values = policy.forward_batched(pf_t, gf_t, mask_t)  # [B, P_max], [B]
        dist = Categorical(logits=logits)
        new_logps = dist.log_prob(actions)
        entropies = dist.entropy()

        ratio = (new_logps - old_logps).exp()
        surr1 = ratio * advs_t
        surr2 = torch.clamp(ratio, 1 - clip, 1 + clip) * advs_t
        policy_loss = -torch.min(surr1, surr2).mean()
        value_loss = F.mse_loss(values, returns_t)
        ent = entropies.mean()
        loss = policy_loss + vc * value_loss - ec * ent

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(policy.parameters(), 0.5)
        optimizer.step()

        with torch.no_grad():
            approx_kl = (old_logps - new_logps).mean().item()
            # Exclude -inf padded positions before computing max magnitude.
            valid_logits = logits.detach().masked_fill(~mask_t, 0.0)
            logit_max = float(valid_logits.abs().max().item())
        last_pl, last_vl, last_ent = policy_loss.item(), value_loss.item(), ent.item()
        last_kl = approx_kl
        last_logit_max = logit_max
    return last_pl, last_vl, last_ent, last_kl, last_logit_max


def greedy_rollout(env: GbrlEnv, policy: PointerPolicy, instance: str, device,
                   max_steps: int, max_wall_s: float = 60.0):
    """Greedy (argmax) rollout. Bails out early if wall-clock exceeds `max_wall_s`
    — eval on cyclic-6+ can hit pathological basis growth, so an open-ended
    greedy eval can take 30+ CPU-min without this cap.

    Returns (steps, basis_size, arith_ops, truncated, wall_s). `truncated` is
    True if we bailed early via max_steps or max_wall_s — i.e. the basis
    isn't a complete Gröbner basis and the reported metrics are partial."""
    import time as _time
    obs = env.reset(instance)
    steps = 0
    t0 = _time.time()
    truncated = False
    while True:
        f = featurize(obs)
        if f is None:
            break
        if _time.time() - t0 > max_wall_s:
            truncated = True
            break
        pair_feats, global_feats = f
        pf = torch.from_numpy(pair_feats).to(device)
        gf = torch.from_numpy(global_feats).to(device)
        with torch.no_grad():
            logits, _ = policy(pf, gf)
            a = int(torch.argmax(logits).item())
        result = env.step(a)
        steps += 1
        obs = result.obs
        if result.done:
            break
        if steps >= max_steps:
            truncated = True
            break
    wall = _time.time() - t0
    basis = obs.get("basis_size", 0)
    arith_ops = obs.get("arith_ops", 0)
    return steps, basis, arith_ops, truncated, wall


def main():
    import time as _time
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", default=os.path.join(os.path.dirname(__file__), "..", "target", "release", "rl_server"))
    # Instances: --instance (singular) and --instances (plural) both accepted.
    ap.add_argument("--instances", nargs="+", default=None)
    ap.add_argument("--instance", default=None, help="alias for --instances with one entry")
    # Stop condition: --updates (preferred) or --episodes.
    ap.add_argument("--updates", type=int, default=None,
                    help="stop after this many PPO updates. Overrides --episodes.")
    ap.add_argument("--episodes", type=int, default=200)
    ap.add_argument("--max-steps", type=int, default=2000)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--gamma", type=float, default=0.99)
    ap.add_argument("--lam", type=float, default=0.95)
    ap.add_argument("--clip", type=float, default=0.2)
    ap.add_argument("--vc", type=float, default=0.5)
    ap.add_argument("--ec", type=float, default=0.01, help="entropy bonus coefficient")
    ap.add_argument("--entropy", type=float, default=None, help="alias for --ec")
    ap.add_argument("--ppo-epochs", type=int, default=4)
    ap.add_argument("--batch-episodes", type=int, default=8,
                    help="number of trajectories collected per PPO update")
    ap.add_argument("--parallel", type=int, default=1)
    ap.add_argument("--step-timeout", type=float, default=0.0)
    ap.add_argument("--eval-max-wall", type=float, default=60.0)
    ap.add_argument("--reward-mode", choices=["step", "arith_ops"], default="step",
                    help="step: -1 per step (default). arith_ops: -Δ(arith_ops) per step, "
                         "which makes PPO directly optimize the cost metric we report.")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--eval-every", type=int, default=1,
                    help="run a greedy (argmax) eval rollout every N updates")
    ap.add_argument("--save", default=None, help="save final policy weights")
    ap.add_argument("--load", default=None, help="load policy weights at start")
    ap.add_argument("--init", default=None, help="alias for --load")
    args = ap.parse_args()

    # Resolve aliases.
    if args.instance is not None:
        args.instances = [args.instance]
    if args.instances is None:
        args.instances = ["cyclic-4"]
    if args.entropy is not None:
        args.ec = args.entropy
    if args.init is not None and args.load is None:
        args.load = args.init

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    device = torch.device(args.device)
    binary_path = os.path.abspath(args.binary)
    if not os.path.exists(binary_path):
        raise SystemExit(f"binary not found: {binary_path}\nbuild first: cargo build --release")

    vec_env = VectorEnv(binary_path, args.parallel)
    eval_env = GbrlEnv(binary_path)  # separate single-env subprocess for greedy eval

    policy = PointerPolicy().to(device)
    if args.load and os.path.exists(args.load):
        policy.load_state_dict(torch.load(args.load, map_location=device))
        print(f"loaded weights from {args.load}", flush=True)
    optimizer = torch.optim.Adam(policy.parameters(), lr=args.lr)

    timeout_str = f"{args.step_timeout}s" if args.step_timeout > 0 else "off"
    target_updates = args.updates if args.updates is not None else 10**9
    target_episodes = args.episodes if args.updates is None else 10**9
    print(f"instances={args.instances}  parallel={args.parallel}  "
          f"batch_ep={args.batch_episodes}  max_steps={args.max_steps}  "
          f"step_timeout={timeout_str}  eval_max_wall={args.eval_max_wall}s  "
          f"lr={args.lr}  ec={args.ec}  clip={args.clip}  ppo_epochs={args.ppo_epochs}  "
          f"reward={args.reward_mode}",
          flush=True)
    if args.updates is not None:
        print(f"stop after {args.updates} updates", flush=True)
    else:
        print(f"stop after {args.episodes} episodes", flush=True)

    history = {inst: deque(maxlen=40) for inst in args.instances}
    updates = 0
    total_episodes = 0
    timeout_s = args.step_timeout if args.step_timeout > 0 else None
    # Best-untruncated tracking: lower arith_ops wins, ties broken by lower steps.
    best_ops = float("inf")
    best_steps = float("inf")
    best_path = (args.save + ".best") if args.save else None
    train_t0 = _time.time()
    try:
        while updates < target_updates and total_episodes < target_episodes:
            upd_t0 = _time.time()
            completed, rmetrics = rollout_async_parallel(
                vec_env, policy, args.instances, args.max_steps,
                args.batch_episodes, device, step_timeout_s=timeout_s,
                reward_mode=args.reward_mode)
            rollout_wall = _time.time() - upd_t0

            batch = []
            batch_steps = []
            for inst, traj, steps, _status in completed:
                total_episodes += 1
                if not traj:
                    continue
                advs, returns = gae(traj, args.gamma, args.lam)
                for k, t in enumerate(traj):
                    t["adv"] = advs[k]
                    t["return"] = returns[k]
                batch.extend(traj)
                batch_steps.append(steps)
                history[inst].append(steps)
            if not batch:
                continue

            advs_arr = np.array([b["adv"] for b in batch], dtype=np.float32)
            advs_arr = (advs_arr - advs_arr.mean()) / (advs_arr.std() + 1e-8)
            for i, b in enumerate(batch):
                b["adv"] = float(advs_arr[i])

            update_t0 = _time.time()
            pl, vl, ent, kl, logit_max = ppo_update_batched(
                policy, optimizer, batch, device,
                args.ppo_epochs, args.clip, args.vc, args.ec)
            update_wall = _time.time() - update_t0
            updates += 1

            B = max(1, args.batch_episodes)
            timeout_rate = rmetrics["n_timeout"] / B
            trunc_rate = rmetrics["n_truncated"] / B
            comp_rate = rmetrics["n_completed"] / B
            all_lens = rmetrics["all_lens"]
            comp_lens = rmetrics["completed_lens"]
            mean_all = float(np.mean(all_lens)) if all_lens else float("nan")
            mean_comp = float(np.mean(comp_lens)) if comp_lens else float("nan")
            upd_wall = _time.time() - upd_t0
            avg = {k: f"{np.mean(v):5.1f}" for k, v in history.items() if v}

            print(
                f"upd={updates:4d} ep={total_episodes:5d} "
                f"len_all={mean_all:6.1f} len_comp={mean_comp:6.1f} "
                f"comp/trunc/to={comp_rate:.2f}/{trunc_rate:.2f}/{timeout_rate:.2f} "
                f"pl={pl:+.4f} vl={vl:7.2f} ent={ent:.2f} kl={kl:+.3f} "
                f"|logit|_max={logit_max:5.2f} "
                f"wall={upd_wall:5.1f}s (roll {rollout_wall:5.1f}s, upd {update_wall:4.1f}s) "
                f"avg={avg}",
                flush=True,
            )

            if updates % args.eval_every == 0:
                for inst in args.instances:
                    steps, basis, ops, trunc, wall = greedy_rollout(
                        eval_env, policy, inst, device,
                        args.max_steps, max_wall_s=args.eval_max_wall)
                    tag = "TRUNC" if trunc else "OK"
                    print(
                        f"  [eval] {tag} {inst}: steps={steps} basis={basis} "
                        f"arith_ops={ops} wall={wall:.1f}s",
                        flush=True,
                    )
                    if not trunc and ops < best_ops or (ops == best_ops and steps < best_steps):
                        best_ops, best_steps = ops, steps
                        if best_path:
                            torch.save(policy.state_dict(), best_path)
                            print(
                                f"  [eval] new best (untruncated): arith_ops={ops} "
                                f"steps={steps} → {best_path}", flush=True,
                            )
    finally:
        train_wall = _time.time() - train_t0
        print(f"\ntotal training wall: {train_wall:.1f}s   "
              f"updates: {updates}   episodes: {total_episodes}", flush=True)
        if math.isfinite(best_ops):
            print(f"best untruncated greedy: arith_ops={best_ops} steps={best_steps}", flush=True)
        else:
            print("best untruncated greedy: NONE (every eval truncated)", flush=True)
        if args.save:
            torch.save(policy.state_dict(), args.save)
            print(f"saved final weights to {args.save}", flush=True)
        vec_env.close()
        eval_env.close()


if __name__ == "__main__":
    main()
