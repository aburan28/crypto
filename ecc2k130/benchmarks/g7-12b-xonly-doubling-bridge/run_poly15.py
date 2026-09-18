#!/usr/bin/env python3
"""Build, validate, time and summarize 15 B/s arithmetic-only and poly12 maps."""
import subprocess, sys
from pathlib import Path

OUT = Path(__file__).resolve().parent
STEPS = [
    "build_poly15.py",
    "validate_poly15.py",
    "screen_poly15.py",
    "summarize_poly15.py",
]


def main():
    for name in STEPS:
        print(f"== {name} ==", flush=True)
        result = subprocess.run([sys.executable, str(OUT / name)], cwd=OUT)
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
