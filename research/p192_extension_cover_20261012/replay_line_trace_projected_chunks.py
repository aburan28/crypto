"""Resume independent per-pivot pure-Python replays of complete chunks."""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("--chunks", type=Path, required=True)
parser.add_argument("--start-pivot", type=int, default=0)
parser.add_argument("--end-pivot", type=int, default=141)
parser.add_argument("--available-only", action="store_true")
parser.add_argument("--resume", action="store_true")
args = parser.parse_args()
assert 0 <= args.start_pivot < args.end_pivot <= 141
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
verifier = HERE / "verify_line_trace_projected_chunk.py"
verifier_hash = sha(verifier)
done = skipped = 0
for i in range(args.start_pivot, args.end_pivot):
    chunk = args.chunks / f"pivot-{i:03d}.json"
    if not chunk.is_file():
        if args.available_only:
            print("TRACE_REPLAY_WAITING_FOR_PIVOT", i, flush=True)
            break
        raise FileNotFoundError(chunk)
    out = args.chunks / f"pivot-{i:03d}-replay.json"
    if out.exists():
        assert args.resume, f"existing replay: {out}"
        prior = json.loads(out.read_text())
        assert prior["status"] == "verified_independent_p192_trace_projected_pivot_replay"
        assert prior["pivot"] == i and prior["chunk_sha256"] == sha(chunk)
        assert prior["source_sha256"] == verifier_hash
        skipped += 1
        continue
    log = args.chunks / f"pivot-{i:03d}-replay.log"
    stderr = args.chunks / f"pivot-{i:03d}-replay.stderr"
    command = [sys.executable, str(verifier), "--chunk", str(chunk),
               "--out", str(out)]
    with log.open("wb") as stdout_file, stderr.open("wb") as stderr_file:
        result = subprocess.run(command, cwd=Path.cwd(), stdout=stdout_file,
                                stderr=stderr_file, check=False)
    assert result.returncode == 0, f"replay failed for pivot {i}: {stderr}"
    assert out.is_file()
    done += 1
    print("TRACE_REPLAY_PIVOT_PASS", i, flush=True)
print("TRACE_REPLAY_DRIVER_PASS", json.dumps({
    "new_replays": done,
    "skipped_verified_replays": skipped,
    "available_only": args.available_only,
}, sort_keys=True), flush=True)
