"""Prepare with builds confined to CPUs 0–1 and measurement CPU 7 available."""
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import tournament

original = tournament.snapshot_build


def confined_build(source, destination):
    affinity = os.sched_getaffinity(0)
    try:
        os.sched_setaffinity(0, {0, 1})
        return original(source, destination)
    finally:
        os.sched_setaffinity(0, affinity)


tournament.snapshot_build = confined_build
tournament.main()
