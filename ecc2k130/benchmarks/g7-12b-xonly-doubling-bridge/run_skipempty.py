#!/usr/bin/env python3
"""Build, validate, time and summarize the empty-queue skip on this G7 GPU."""
import subprocess, sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
STEPS = [
    "build_skipempty.py",
    "validate_skipempty.py",
    "screen_skipempty.py",
    "screen_mod72_dp34.py",
    "summarize_skipempty.py",
]


def main():
    for name in STEPS:
        print(f"== {name} ==", flush=True)
        result = subprocess.run([sys.executable, str(OUT / name)], cwd=OUT)
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
