#!/bin/sh
# Stage gpu/ecc2k under this directory so ecc2k130/modal_job.py (which
# uploads only the ecc2k130 tree) can build and run it:
#   cd ecc2k130 && benchmarks/ecc2k-clmad/prepare.sh && \
#   modal run --detach modal_job.py --job benchmarks/ecc2k-clmad/gpujob.sh --out /tmp/ecc2k-clmad
# tree/ is a copy, ignored by git; the sources of record are gpu/ecc2k.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
rsync -a --delete --exclude 'test_*' --exclude 'bench' --exclude '*.bin' --exclude '*.log' \
  --exclude 'curve_*.h' --exclude 'vec_*.h' --exclude '__pycache__' \
  "$here/../../../gpu/ecc2k/" "$here/tree/"
