"""Compare a fresh bounded run with frozen evidence, ignoring only elapsed time."""
import gzip
import json
import sys
from pathlib import Path


def read(path):
    data = json.loads(gzip.decompress(path.read_bytes()))
    del data['elapsed_seconds']
    return data


if __name__ == '__main__':
    frozen = Path(__file__).resolve().parent/'results'
    fresh = Path(sys.argv[1])
    if read(frozen/'raw.json.gz') != read(fresh/'raw.json.gz'):
        raise SystemExit('raw evidence differs')
    if json.loads((frozen/'summary.json').read_text()) != json.loads((fresh/'summary.json').read_text()):
        raise SystemExit('summary differs')
    print('Exact replay: all inputs, certificates, roots, traces and counters match.')
