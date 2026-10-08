"""
End-to-end demo: generate a prime-order curve over a `bits`-bit prime, pick a
random secret k, Q = k*P, recover k with the batched rho harness.

  cd tpu && python3 -m rho.demo [bits] [numWalks] [dpBits]

Correctness demo only. Steps/s printed by the solver are CPU-backend numbers
and are not a performance claim for any device.
"""

import random
import sys
import time

from rho import curve as curveLib
from rho.solve import RhoSolver


def main():
    bits = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    numWalks = int(sys.argv[2]) if len(sys.argv) > 2 else 2048
    dpBits = int(sys.argv[3]) if len(sys.argv) > 3 else 8

    rng = random.Random(2026)
    curve = curveLib.generatePrimeOrderCurve(bits, rng)
    secret = rng.randrange(1, curve.order)
    Q = curve.mul(secret, curve.generator)
    print(f"p = {curve.p}\nn = {curve.order}\nsecret k = {secret}")

    solver = RhoSolver(curve, Q, numWalks=numWalks, stepsPerChunk=32, tableSize=32, dpBits=dpBits, seed=1)
    startTime = time.time()
    k = solver.solve()
    print(f"recovered k = {k}  ({'correct' if k == secret else 'WRONG'}), "
          f"{time.time() - startTime:.1f}s wall, {solver.badCount} bad-denominator restarts")


if __name__ == "__main__":
    main()
