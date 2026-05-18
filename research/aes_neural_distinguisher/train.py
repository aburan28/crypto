"""Train a Gohr-style neural distinguisher on round-reduced AES-128.

Example:
    python train.py --rounds 3 --delta 80000000000000000000000000000000 \\
        --epochs 20 --batch 5000 --batches-per-epoch 1000

The default delta activates a single byte in column 0, which is a reasonable
seed for 3-round AES — past round 1 SubBytes the active byte spreads to a
full column via MixColumns, and after round 2 to the whole state. Replace
with a delta from a MILP-found truncated differential for harder round counts.
"""

from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from data import PairStream, fixed_eval_set
from model import build_model, param_count


def matthews_corrcoef(pred: torch.Tensor, y: torch.Tensor) -> float:
    p = (pred >= 0.5).float()
    tp = ((p == 1) & (y == 1)).sum().item()
    tn = ((p == 0) & (y == 0)).sum().item()
    fp = ((p == 1) & (y == 0)).sum().item()
    fn = ((p == 0) & (y == 1)).sum().item()
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return 0.0 if denom == 0 else (tp * tn - fp * fn) / denom


@torch.no_grad()
def evaluate(model: nn.Module, X: torch.Tensor, y: torch.Tensor, device: str, batch: int = 4096):
    model.eval()
    probs = []
    for i in range(0, X.shape[0], batch):
        xb = X[i : i + batch].to(device)
        probs.append(torch.sigmoid(model(xb)).cpu())
    p = torch.cat(probs)
    acc = ((p >= 0.5).float() == y).float().mean().item()
    mcc = matthews_corrcoef(p, y)
    return acc, mcc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--delta",  type=str, default="80000000000000000000000000000000",
                    help="32-hex-char input difference (default: 0x80 in byte 0).")
    ap.add_argument("--encoding", choices=["byte", "bits"], default="bits",
                    help="Input encoding. 'bits' is the Gohr-style 128-bit representation.")
    ap.add_argument("--task", choices=["random", "wrong_delta"], default="random",
                    help="Negative-class generator. 'wrong_delta' = Gohr's real-difference setup.")
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--batch",  type=int, default=5000)
    ap.add_argument("--batches-per-epoch", type=int, default=1000)
    ap.add_argument("--eval-size", type=int, default=100_000)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--weight-decay", type=float, default=1e-5)
    ap.add_argument("--scheduler", choices=["onecycle", "constant"], default="onecycle")
    ap.add_argument("--channels", type=int, default=64)
    ap.add_argument("--blocks", type=int, default=8)
    ap.add_argument("--num-workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", type=str, default=None)
    ap.add_argument("--save", type=str, default=None)
    ap.add_argument("--log-file", type=str, default=None,
                    help="Optional line-oriented training log for long background runs.")
    ap.add_argument("--init-from", type=str, default=None,
                    help="Path to a checkpoint to warm-start model weights from.")
    args = ap.parse_args()

    log_fh = None
    if args.log_file:
        log_path = Path(args.log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_fh = log_path.open("a", buffering=1)

    def emit(msg: str) -> None:
        print(msg, flush=True)
        if log_fh is not None:
            log_fh.write(msg + "\n")

    device = args.device or ("cuda" if torch.cuda.is_available()
                             else "mps" if torch.backends.mps.is_available()
                             else "cpu")
    torch.manual_seed(args.seed)

    model = build_model(args.encoding, channels=args.channels, blocks=args.blocks).to(device)
    if args.init_from:
        ckpt = torch.load(args.init_from, map_location=device, weights_only=False)
        missing, unexpected = model.load_state_dict(ckpt["model"], strict=False)
        emit(f"warm-start from {args.init_from}  missing={len(missing)}  unexpected={len(unexpected)}")
    emit(f"device={device}  params={param_count(model):,}  rounds={args.rounds}  "
         f"delta={args.delta}  encoding={args.encoding}  task={args.task}")

    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    total_steps = args.epochs * args.batches_per_epoch
    if args.scheduler == "onecycle":
        sched = torch.optim.lr_scheduler.OneCycleLR(
            opt, max_lr=args.lr, total_steps=total_steps, pct_start=0.1, anneal_strategy="cos"
        )
    else:
        sched = torch.optim.lr_scheduler.ConstantLR(opt, factor=1.0, total_iters=1)
    loss_fn = nn.BCEWithLogitsLoss()

    eval_X, eval_y = fixed_eval_set(args.rounds, args.delta, args.eval_size,
                                    encoding=args.encoding, task=args.task)
    emit(f"eval set: {args.eval_size} pairs, fixed seed")

    train_ds = PairStream(args.rounds, args.delta, args.batch, args.batches_per_epoch,
                          encoding=args.encoding, task=args.task, seed=args.seed + 1)
    loader = DataLoader(train_ds, batch_size=None, num_workers=args.num_workers,
                        persistent_workers=args.num_workers > 0)

    best_acc = 0.0
    for epoch in range(args.epochs):
        model.train()
        t0 = time.time()
        running, n = 0.0, 0
        for X, y in loader:
            X, y = X.to(device, non_blocking=True), y.to(device, non_blocking=True)
            opt.zero_grad(set_to_none=True)
            logits = model(X)
            loss = loss_fn(logits, y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            sched.step()
            running += loss.item() * X.shape[0]
            n += X.shape[0]
        acc, mcc = evaluate(model, eval_X, eval_y, device)
        dt = time.time() - t0
        emit(f"epoch {epoch+1:3d}/{args.epochs}  loss={running/n:.4f}  acc={acc:.4f}  mcc={mcc:.4f}  ({dt:.1f}s)")
        if acc > best_acc:
            best_acc = acc
            if args.save:
                torch.save({"model": model.state_dict(), "args": vars(args), "acc": acc}, args.save)

    emit(f"best eval acc: {best_acc:.4f}")
    if log_fh is not None:
        log_fh.close()


if __name__ == "__main__":
    main()
