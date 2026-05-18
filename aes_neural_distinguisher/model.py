"""ResNet distinguishers for round-reduced AES.

Distinguisher2D : input (B, 3, 4, 4), 2D convs over the AES state-as-image.
Distinguisher1D : input (B, 3, 128),  1D convs over the bit string (Gohr-style).

Both end in a small MLP head producing a single logit.
"""

from __future__ import annotations

import torch
import torch.nn as nn


class ResBlock2D(nn.Module):
    def __init__(self, ch: int):
        super().__init__()
        self.conv1 = nn.Conv2d(ch, ch, 3, padding=1, bias=False)
        self.bn1   = nn.BatchNorm2d(ch)
        self.conv2 = nn.Conv2d(ch, ch, 3, padding=1, bias=False)
        self.bn2   = nn.BatchNorm2d(ch)
        self.act   = nn.ReLU(inplace=True)

    def forward(self, x):
        r = self.act(self.bn1(self.conv1(x)))
        r = self.bn2(self.conv2(r))
        return self.act(x + r)


class ResBlock1D(nn.Module):
    def __init__(self, ch: int, kernel: int = 3):
        super().__init__()
        pad = kernel // 2
        self.conv1 = nn.Conv1d(ch, ch, kernel, padding=pad, bias=False)
        self.bn1   = nn.BatchNorm1d(ch)
        self.conv2 = nn.Conv1d(ch, ch, kernel, padding=pad, bias=False)
        self.bn2   = nn.BatchNorm1d(ch)
        self.act   = nn.ReLU(inplace=True)

    def forward(self, x):
        r = self.act(self.bn1(self.conv1(x)))
        r = self.bn2(self.conv2(r))
        return self.act(x + r)


class Distinguisher2D(nn.Module):
    def __init__(self, channels: int = 64, blocks: int = 8, mlp_hidden: int = 256):
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv2d(3, channels, 3, padding=1, bias=False),
            nn.BatchNorm2d(channels),
            nn.ReLU(inplace=True),
        )
        self.trunk = nn.Sequential(*[ResBlock2D(channels) for _ in range(blocks)])
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.Linear(channels * 4 * 4, mlp_hidden),
            nn.ReLU(inplace=True),
            nn.Linear(mlp_hidden, mlp_hidden),
            nn.ReLU(inplace=True),
            nn.Linear(mlp_hidden, 1),
        )

    def forward(self, x):
        return self.head(self.trunk(self.stem(x))).squeeze(-1)


class Distinguisher1D(nn.Module):
    def __init__(self, channels: int = 64, blocks: int = 8, mlp_hidden: int = 256, width: int = 128):
        super().__init__()
        self.stem = nn.Sequential(
            nn.Conv1d(3, channels, 1, bias=False),
            nn.BatchNorm1d(channels),
            nn.ReLU(inplace=True),
        )
        self.trunk = nn.Sequential(*[ResBlock1D(channels) for _ in range(blocks)])
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.Linear(channels * width, mlp_hidden),
            nn.ReLU(inplace=True),
            nn.Linear(mlp_hidden, mlp_hidden),
            nn.ReLU(inplace=True),
            nn.Linear(mlp_hidden, 1),
        )

    def forward(self, x):
        return self.head(self.trunk(self.stem(x))).squeeze(-1)


def build_model(encoding: str, channels: int = 64, blocks: int = 8, mlp_hidden: int = 256) -> nn.Module:
    if encoding == "byte":
        return Distinguisher2D(channels, blocks, mlp_hidden)
    if encoding == "bits":
        return Distinguisher1D(channels, blocks, mlp_hidden)
    raise ValueError(f"unknown encoding: {encoding!r}")


def param_count(m: nn.Module) -> int:
    return sum(p.numel() for p in m.parameters() if p.requires_grad)
