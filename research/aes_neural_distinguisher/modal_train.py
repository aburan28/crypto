"""Run the AES distinguisher trainer on a Modal CUDA worker.

Examples:
    modal run modal_train.py --epochs 200 --batches-per-epoch 50 --gpu A10G
    modal volume get aes-neural-distinguisher-runs r3_single_byte_cuda.pt .
    modal volume get aes-neural-distinguisher-runs r3_single_byte_cuda.log .
"""

from __future__ import annotations

import shlex
import sys
from pathlib import Path

import modal


APP_NAME = "aes-neural-distinguisher"
VOLUME_NAME = "aes-neural-distinguisher-runs"

app = modal.App(APP_NAME)
runs = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

image = (
    modal.Image.from_registry("pytorch/pytorch:2.5.1-cuda12.4-cudnn9-runtime")
    .pip_install("numpy")
    .add_local_python_source("aes", "data", "model", "train")
)


@app.function(
    image=image,
    gpu="A10G",
    timeout=24 * 60 * 60,
    volumes={"/runs": runs},
)
def train_remote(
    rounds: int,
    delta: str,
    encoding: str,
    task: str,
    epochs: int,
    batch: int,
    batches_per_epoch: int,
    eval_size: int,
    lr: float,
    weight_decay: float,
    scheduler: str,
    channels: int,
    blocks: int,
    seed: int,
    run_name: str,
    init_from: str | None,
):
    save_path = f"/runs/{run_name}.pt"
    log_path = f"/runs/{run_name}.log"
    init_path = None
    if init_from:
        init_path = init_from if init_from.startswith("/") else f"/runs/{init_from}"
    argv = [
        "train.py",
        "--rounds", str(rounds),
        "--delta", delta,
        "--encoding", encoding,
        "--task", task,
        "--epochs", str(epochs),
        "--batch", str(batch),
        "--batches-per-epoch", str(batches_per_epoch),
        "--eval-size", str(eval_size),
        "--lr", str(lr),
        "--weight-decay", str(weight_decay),
        "--scheduler", scheduler,
        "--channels", str(channels),
        "--blocks", str(blocks),
        "--num-workers", "2",
        "--seed", str(seed),
        "--device", "cuda",
        "--save", save_path,
        "--log-file", log_path,
    ]
    if init_path:
        argv.extend(["--init-from", init_path])
    print("remote command:", " ".join(shlex.quote(x) for x in argv), flush=True)
    sys.argv = argv
    import train

    try:
        train.main()
    finally:
        runs.commit()
    return {"checkpoint": save_path, "log": log_path}


@app.local_entrypoint()
def main(
    rounds: int = 3,
    delta: str = "80000000000000000000000000000000",
    encoding: str = "bits",
    task: str = "wrong_delta",
    epochs: int = 200,
    batch: int = 5000,
    batches_per_epoch: int = 50,
    eval_size: int = 100000,
    lr: float = 3e-4,
    weight_decay: float = 1e-5,
    scheduler: str = "constant",
    channels: int = 48,
    blocks: int = 6,
    seed: int = 42,
    run_name: str = "r3_single_byte_cuda",
    init_from: str | None = None,
    spawn: bool = False,
):
    kwargs = dict(
        rounds=rounds,
        delta=delta,
        encoding=encoding,
        task=task,
        epochs=epochs,
        batch=batch,
        batches_per_epoch=batches_per_epoch,
        eval_size=eval_size,
        lr=lr,
        weight_decay=weight_decay,
        scheduler=scheduler,
        channels=channels,
        blocks=blocks,
        seed=seed,
        run_name=run_name,
        init_from=init_from,
    )
    if spawn:
        call = train_remote.spawn(**kwargs)
        print(f"spawned background call: {call.object_id}")
        print(f"checkpoint: /runs/{run_name}.pt")
        print(f"log       : /runs/{run_name}.log")
        print(f"download : modal volume get {VOLUME_NAME} {run_name}.pt .")
        print(f"download : modal volume get {VOLUME_NAME} {run_name}.log .")
        return

    result = train_remote.remote(**kwargs)
    print(f"checkpoint: {result['checkpoint']}")
    print(f"log       : {result['log']}")
    print(f"download : modal volume get {VOLUME_NAME} {Path(result['checkpoint']).name} .")
    print(f"download : modal volume get {VOLUME_NAME} {Path(result['log']).name} .")
