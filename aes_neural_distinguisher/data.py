"""Streaming pair generator for Gohr-style training.

Label 1 = "real":   C0 = E_K(P), C1 = E_K(P xor delta_in), random K, random P.
Label 0 = "random": C0, C1 independent uniform 128-bit.

Two encodings are supported:
  encoding="byte" -> shape (B, 3, 4, 4), values in [0, 1]   (state-as-image)
  encoding="bits" -> shape (B, 3, 128),  values in {0, 1}   (Gohr-style bits)

Channels are always (C0, C1, C0 ^ C1).
"""

from __future__ import annotations

import numpy as np
import torch
from torch.utils.data import IterableDataset

from aes import encrypt, bytes_to_state


def parse_delta(spec: str) -> np.ndarray:
    spec = spec.replace("_", "").replace(" ", "")
    if len(spec) != 32:
        raise ValueError(f"delta must be 32 hex chars, got {len(spec)}")
    return np.frombuffer(bytes.fromhex(spec), dtype=np.uint8)


def _encode(c0: np.ndarray, c1: np.ndarray, encoding: str) -> np.ndarray:
    if encoding == "byte":
        s0 = bytes_to_state(c0).astype(np.float32) / 255.0
        s1 = bytes_to_state(c1).astype(np.float32) / 255.0
        sx = bytes_to_state(c0 ^ c1).astype(np.float32) / 255.0
        return np.stack([s0, s1, sx], axis=1)                                # (B, 3, 4, 4)
    if encoding == "bits":
        b0 = np.unpackbits(c0, axis=-1, bitorder="big").astype(np.float32)   # (B, 128)
        b1 = np.unpackbits(c1, axis=-1, bitorder="big").astype(np.float32)
        bx = np.unpackbits(c0 ^ c1, axis=-1, bitorder="big").astype(np.float32)
        return np.stack([b0, b1, bx], axis=1)                                # (B, 3, 128)
    raise ValueError(f"unknown encoding: {encoding!r}")


def make_pairs(batch: int, rounds: int, delta: np.ndarray, rng: np.random.Generator,
               encoding: str = "byte", task: str = "random"):
    """task='random'      : label-0 = uniform random ciphertext pair.
       task='wrong_delta' : label-0 = real ciphertext pair from a uniformly random
                            non-zero input difference (Gohr's 'real difference' setup).
                            Forces the model to learn the differential structure of
                            `delta` specifically, not generic real-vs-random features.
    """
    half = batch // 2
    keys = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
    p0 = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
    p1 = p0 ^ delta
    c0_real = encrypt(p0, keys, rounds)
    c1_real = encrypt(p1, keys, rounds)

    if task == "random":
        c0_neg = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
        c1_neg = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
    elif task == "wrong_delta":
        wrong = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
        # Collision with target delta has probability 2^-128 — skip the check.
        # Force non-zero by setting a random byte if a row is all-zero (negligible chance).
        zero_rows = np.where(np.all(wrong == 0, axis=1))[0]
        if zero_rows.size:
            wrong[zero_rows, 0] = rng.integers(1, 256, size=zero_rows.size, dtype=np.uint8)
        keys_w = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
        p0_w = rng.integers(0, 256, size=(half, 16), dtype=np.uint8)
        p1_w = p0_w ^ wrong
        c0_neg = encrypt(p0_w, keys_w, rounds)
        c1_neg = encrypt(p1_w, keys_w, rounds)
    else:
        raise ValueError(f"unknown task: {task!r}")

    c0 = np.concatenate([c0_real, c0_neg], axis=0)
    c1 = np.concatenate([c1_real, c1_neg], axis=0)
    y  = np.concatenate([np.ones(half, dtype=np.float32),
                         np.zeros(half, dtype=np.float32)])
    perm = rng.permutation(batch)
    c0, c1, y = c0[perm], c1[perm], y[perm]
    X = _encode(c0, c1, encoding)
    return X, y


class PairStream(IterableDataset):
    def __init__(self, rounds: int, delta_hex: str, batch: int, batches_per_epoch: int,
                 encoding: str = "byte", task: str = "random", seed: int = 0):
        self.rounds = rounds
        self.delta = parse_delta(delta_hex)
        self.batch = batch
        self.batches_per_epoch = batches_per_epoch
        self.encoding = encoding
        self.task = task
        self.seed = seed

    def __iter__(self):
        info = torch.utils.data.get_worker_info()
        worker_id = 0 if info is None else info.id
        rng = np.random.default_rng(self.seed + worker_id * 99991)
        for _ in range(self.batches_per_epoch):
            X, y = make_pairs(self.batch, self.rounds, self.delta, rng, self.encoding, self.task)
            yield torch.from_numpy(X), torch.from_numpy(y)


def fixed_eval_set(rounds: int, delta_hex: str, n: int, encoding: str = "byte",
                   task: str = "random", seed: int = 12345):
    rng = np.random.default_rng(seed)
    delta = parse_delta(delta_hex)
    X, y = make_pairs(n, rounds, delta, rng, encoding, task)
    return torch.from_numpy(X), torch.from_numpy(y)
