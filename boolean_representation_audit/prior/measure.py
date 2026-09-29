"""Declared synthetic controls. No target data or cryptographic integration."""
import hashlib
import json
import platform
import random
import statistics
import time
import tracemalloc
from pathlib import Path

from boolean_closure import compute, equal_spans, solutions, verify


def cases():
    for n in (6, 8):
        yield f"monomial-{n}", n, [[3]]
        yield f"block-linear-{n}", n, [[1 << i, 1 << (i + 1)]
                                     for i in range(0, n, 2)]
        rng = random.Random(131 + n)
        pool = [m for m in range(1, 1 << n) if m.bit_count() <= 2]
        point = rng.randrange(1 << n)
        planted = []
        for _ in range(n // 2):
            row = rng.sample(pool, 5)
            if sum((m & point) == m for m in row) % 2:
                row.append(0)
            planted.append(row)
        yield f"planted-quadratic-{n}", n, planted
        yield f"dependent-contradictory-{n}", n, [[1], [1], [1, 0]]


def main():
    records = []
    for name, n, generators in cases():
        entry = {"name": name, "n": n, "generators": generators,
                 "solution_count": len(solutions(n, generators)), "algorithms": {}}
        complete_results = {}
        for strategy in ("exhaustive", "frontier"):
            times = []
            for _ in range(3):
                start = time.perf_counter()
                result = compute(n, generators, strategy=strategy)
                times.append(time.perf_counter() - start)
                verify(n, generators, result)
            tracemalloc.start()
            traced = compute(n, generators, strategy=strategy)
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            assert result == traced
            assert result["stats"]["rank"] == (1 << n) - entry["solution_count"]
            entry["algorithms"][strategy] = {
                "wall_seconds": times, "median_seconds": statistics.median(times),
                "peak_python_bytes": peak, "stats": result["stats"],
                "verification": verify(n, generators, result)}
            complete_results[strategy] = result
        assert equal_spans(*complete_results.values())
        entry["equal_spans"] = True
        records.append(entry)
        print(name, {s: {"seconds": round(v["median_seconds"], 5),
                        "rows": v["stats"]["submitted"],
                        "peak_bytes": v["peak_python_bytes"]}
                     for s, v in entry["algorithms"].items()}, flush=True)
    payload = {"protocol_sha256": hashlib.sha256(Path("PROTOCOL.md").read_bytes()).hexdigest(),
               "python": platform.python_version(), "platform": platform.platform(),
               "records": records}
    Path("results.json").write_text(json.dumps(payload, indent=2) + "\n")
    example = compute(3, [[1, 0]])
    Path("example_certificate.json").write_text(json.dumps({
        "n": 3, "generators": [[1, 0]], "result": example,
        "verification": verify(3, [[1, 0]], example)}, indent=2) + "\n")


if __name__ == "__main__":
    main()
