"""ECC2K-130 Pollard-rho walk step as a KernelBench problem for KernelAgent.

One step of the canonical ECC2K-130 iteration for a batch of affine points:

    j  = 3 + ((HW(x) >> 1) & 7)          HW = weight of the normal-basis coordinates
    P' = P + sigma^j(P)                  sigma = Frobenius, a coordinate permutation

over GF(2^131) in the permuted type-II optimal normal basis of
include/ref.h (ScalarOnb, NRING = 263), with the device's special-case-free
affine addition (ScalarOnb::addPtRaw):

    d  = x + x2,  lam = (y + y2) / d,  x3 = lam^2 + lam + d,  y3 = lam (x + x3) + x3 + y

Inversion is Itoh-Tsujii exactly as ScalarOnb::inv, so inv(0) = 0 and the
result is defined for every input.  An element is 131 bits in three int64
limbs, bit i of the normal-basis vector at limb i // 64, bit i % 64 (the
device layout); limb 2 holds bits 128..130 and is zero above them.

forward(x, y) takes two (N, 3) int64 tensors and returns an (N, 2, 3) int64
tensor [x3, y3].  This reference is written for clarity, not speed: the
product is the cyclic convolution mod 2 over Z/263, done here with a float64
FFT.  check_reference.py checks it bit for bit against include/ref.h.
"""
import torch
import torch.nn as nn

M = 131
NRING = 2 * M + 1  # 263
LIMBS = 3


def _sigma_tables() -> torch.Tensor:
    """table[k, dst] = src with sigma^k(a)[dst] = a[src], as ScalarOnb::sigma."""
    table = torch.empty(M, M, dtype=torch.long)
    e = 1
    for k in range(M):
        for i in range(1, M + 1):
            t = i * e % NRING
            if t > M:
                t = NRING - t
            table[k, t - 1] = i - 1
        e = 2 * e % NRING
    return table


class Model(nn.Module):
    def __init__(self):
        super().__init__()
        self.register_buffer("sigma_table", _sigma_tables(), persistent=False)
        self.register_buffer("shifts", torch.arange(64, dtype=torch.long), persistent=False)

    # ---- representation ------------------------------------------------
    def unpack(self, a: torch.Tensor) -> torch.Tensor:
        bits = (a.unsqueeze(-1) >> self.shifts) & 1          # (N, 3, 64)
        return bits.reshape(a.shape[0], LIMBS * 64)[:, :M]   # (N, 131)

    def pack(self, bits: torch.Tensor) -> torch.Tensor:
        n = bits.shape[0]
        full = torch.zeros(n, LIMBS * 64, dtype=torch.long, device=bits.device)
        full[:, :M] = bits
        return (full.reshape(n, LIMBS, 64) << self.shifts).sum(-1)  # disjoint bits: sum == or

    # ---- field, on (N, 131) bit tensors --------------------------------
    def sigma(self, a: torch.Tensor, k) -> torch.Tensor:
        if isinstance(k, int):
            return a[:, self.sigma_table[k % M]]
        return a.gather(1, self.sigma_table[k % M])

    def mul(self, a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
        n = a.shape[0]
        lift_a = torch.zeros(n, NRING, dtype=torch.float64, device=a.device)
        lift_b = torch.zeros_like(lift_a)
        lift_a[:, 1:M + 1] = a.double()
        lift_a[:, M + 1:] = a.flip(1).double()               # index NRING - i <- bit i - 1
        lift_b[:, 1:M + 1] = b.double()
        lift_b[:, M + 1:] = b.flip(1).double()
        conv = torch.fft.irfft(torch.fft.rfft(lift_a, dim=1) * torch.fft.rfft(lift_b, dim=1),
                               n=NRING, dim=1)
        acc = torch.round(conv).long() & 1
        acc = acc ^ acc[:, :1]                                # normalise: add all-ones if bit 0
        return acc[:, 1:M + 1]

    def inv(self, a: torch.Tensor) -> torch.Tensor:
        acc, k, e = a, 1, M - 1
        hb = e.bit_length() - 1
        for b in range(hb - 1, -1, -1):
            acc = self.mul(self.sigma(acc, k), acc)
            k <<= 1
            if (e >> b) & 1:
                acc = self.mul(self.sigma(acc, 1), a)
                k += 1
        return self.sigma(acc, 1)

    # ---- the walk step -------------------------------------------------
    def forward(self, x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
        xb, yb = self.unpack(x), self.unpack(y)
        hw = xb.sum(1)
        j = 3 + ((hw >> 1) & 7)
        x2, y2 = self.sigma(xb, j), self.sigma(yb, j)
        d = xb ^ x2
        lam = self.mul(yb ^ y2, self.inv(d))
        x3 = self.sigma(lam, 1) ^ lam ^ d
        y3 = self.mul(lam, xb ^ x3) ^ x3 ^ yb
        return torch.stack([self.pack(x3), self.pack(y3)], dim=1)


BATCH = 1 << 16


def _random_elements(n: int, generator=None) -> torch.Tensor:
    a = torch.randint(-(1 << 63), (1 << 63) - 1, (n, LIMBS), dtype=torch.long, generator=generator)
    a[:, 2] &= 7
    return a


def get_inputs():
    return [_random_elements(BATCH), _random_elements(BATCH)]


def get_init_inputs():
    return []
