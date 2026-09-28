#!/usr/bin/env python3
"""A cairn orbit-piecework job document for one of this client's curves.

`build/witness` emits the claim artifact a cairn node pays for, and the only
way to know the artifact is right is to hand it to cairn's own checker.  That
checker reads a *job document*, which pins the field, the curve, the walk and
the normal basis the orbit is named in.  cairn ships one for ECC2K-130; this
writes the equivalent for the small curves, where a trail is short enough to
walk in the time a test may take, so the whole protocol can be exercised end
to end rather than argued about.

Every constant comes out of the generated header, so a job written here is the
curve the client actually walks:

    python3 cairn_job.py --curve 41 --instance 2 --dp-weight 12 > job.json
    python3 <cairn>/examples/certicom-ecdlp/tools/orbit_dp.py verify \\
        --job job.json claims.json

The normal element is this basis's first conjugate, which is a legitimate
choice for any job -- cairn's own note says any normal basis makes sigma a
rotation -- and is the one `build/witness --nb-generator` takes.  The point of
the exercise is that the two implementations agree on the *name*, not that
they happen to share a generator.
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def header(curve):
    for name in ("eccF%d" % curve,):
        path = os.path.join(HERE, "generated", "%s.h" % name)
        if os.path.exists(path):
            return name, open(path).read()
    sys.exit("no generated header for curve %d; run `make generate`" % curve)


def scalar(src, name, cast=int):
    m = re.search(r"static const (?:int|char \*) ?%s = ([^;]+);" % name, src)
    if not m:
        sys.exit("header has no %s" % name)
    return cast(m.group(1).strip().strip('"'))


def vec(src, name):
    m = re.search(r"%s\[3\] = \{([^}]*)\}" % name, src)
    if not m:
        sys.exit("header has no %s" % name)
    limbs = [int(x.strip().rstrip("ull"), 16) for x in m.group(1).split(",") if x.strip()]
    return limbs[0] | (limbs[1] << 64) | (limbs[2] << 128)


def table(src, name, m):
    mt = re.search(r"%s\[%d\]\[3\] = \{(.*?)\n\};" % (name, m), src, re.S)
    if not mt:
        sys.exit("header has no %s" % name)
    out = []
    for row in re.findall(r"\{([^}]*)\}", mt.group(1)):
        limbs = [int(x.strip().rstrip("ull"), 16) for x in row.split(",") if x.strip()]
        out.append(limbs[0] | (limbs[1] << 64) | (limbs[2] << 128))
    return out


def instance_row(src, name, index, m):
    mt = re.search(r"%s\[\d+\]\[3\] = \{(.*?)\n\};" % name, src, re.S)
    if not mt:
        return None
    rows = re.findall(r"\{([^}]*)\}", mt.group(1))
    if index >= len(rows):
        sys.exit("curve has %d planted instances" % len(rows))
    limbs = [int(x.strip().rstrip("ull"), 16) for x in rows[index].split(",") if x.strip()]
    return limbs[0] | (limbs[1] << 64) | (limbs[2] << 128)


def to_poly(coords, gamma_to_pb, m):
    """ONB coordinates -> the polynomial basis the job names."""
    acc = 0
    for i in range(m):
        if (coords >> i) & 1:
            acc ^= gamma_to_pb[i]
    return acc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curve", type=int, default=41)
    ap.add_argument("--instance", type=int, default=-1)
    ap.add_argument("--dp-weight", type=int, default=None)
    ap.add_argument("--max-batch", type=int, default=64)
    ap.add_argument("--max-steps", type=int, default=1 << 20)
    ap.add_argument("--name", default=None)
    ap.add_argument("--nb-generator", default=None,
                    help="normal element as polynomial-basis hex; must be a conjugate "
                         "of this basis. Defaults to its first, which makes the "
                         "coordinate map trivial -- pass a real job's generator to "
                         "exercise the permutation.")
    args = ap.parse_args()

    ns, src = header(args.curve)
    m = scalar(src, "M")
    taps = re.search(r"PB_TAPS\[3\] = \{([^}]*)\}", src)
    if not taps:
        sys.exit("header has no PB_TAPS")
    tapv = [int(x.strip()) for x in taps.group(1).split(",")]
    poly = (1 << m) | 1
    for t in tapv:
        if t >= 0:
            poly |= 1 << t

    ell = scalar(src, "ELL_DEC", str)
    s = scalar(src, "S_DEC", str)
    dpw = args.dp_weight if args.dp_weight is not None else scalar(src, "DP_WEIGHT")
    gamma_to_pb = table(src, "GAMMA_TO_PB", m)

    if args.instance >= 0:
        px = instance_row(src, "INSTANCE_PX", args.instance, m)
        py = instance_row(src, "INSTANCE_PY", args.instance, m)
        qx = instance_row(src, "INSTANCE_QX", args.instance, m)
        qy = instance_row(src, "INSTANCE_QY", args.instance, m)
        if px is None:
            sys.exit("curve %d has no planted instances" % args.curve)
    else:
        px, py, qx, qy = (vec(src, "PX"), vec(src, "PY"), vec(src, "QX"), vec(src, "QY"))

    if args.nb_generator is not None and int(args.nb_generator, 16) not in gamma_to_pb:
        sys.exit("nb_generator is not a conjugate of this basis: no coordinate "
                 "permutation exists, so witness would emit names nobody can check")

    job = {
        "a": "0",
        "b": "1",
        "cofactor": 4,
        "dp_max_weight": dpw,
        "family": "koblitz-frobenius-rho",
        "frobenius_eigenvalue": "%x" % int(s),
        "generator": {"x": "%x" % to_poly(px, gamma_to_pb, m),
                      "y": "%x" % to_poly(py, gamma_to_pb, m)},
        "j_base": 3,
        "j_count": 8,
        "m": m,
        "max_batch": args.max_batch,
        "max_steps_per_walker": args.max_steps,
        "name": args.name or ("ecc2k-%d client curve" % m),
        "nb_generator": args.nb_generator or ("%x" % gamma_to_pb[0]),
        "order": "%x" % int(ell),
        "poly": "%x" % poly,
        "seed": 0,
        "start_terms": 128,
        "target": {"x": "%x" % to_poly(qx, gamma_to_pb, m),
                   "y": "%x" % to_poly(qy, gamma_to_pb, m)},
        "trail_bits": 16,
        "unit_size": 256,
        "units": 1 << 48,
        "version": 2,
        "witness": "j-counts",
    }
    json.dump(job, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
