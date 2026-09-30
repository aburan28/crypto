"""Use the repository isolation CLI and unwind its affinity restoration on TERM."""
from pathlib import Path
import signal
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[5]/'tools'))
import isolated_bench


def interrupted(signum,frame):
    raise SystemExit(128+signum)


if __name__=='__main__':
    signal.signal(signal.SIGTERM,interrupted)
    raise SystemExit(isolated_bench.main())
