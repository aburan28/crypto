#!/usr/bin/env python3
"""End-to-end test for the orbit witness generator.

Builds a corpus with the client, turns it into claim artifacts with
`build/witness`, and holds the result to the rules a cairn node applies before
it pays.

The client writes a v2 corpus when it is built with the witness counters, so
this covers both ways of getting the counts: straight from a v2 record, and by
replaying a v1 record the way every pre-counter corpus has to be read.  The v1
corpus is derived from the v2 one here, so the two runs describe the same
trails and their `j` vectors must agree exactly -- which is the check that the
counters carried on the device are the same numbers a replay recovers.

The rules:

  * the batch is an object with the single key `dps`, at most `max_batch` long;
  * every element carries exactly `x`, `seed`, `j`;
  * `x` and `seed` are lowercase hex with no prefix and no leading zero;
  * `x` fits `m` normal-basis coordinates and its weight is distinguished;
  * `x` is the *least* rotation of its orbit -- rechecked here from the emitted
    value, so a canonicalisation bug cannot hide behind the tool that made it;
  * `j` is eight non-negative counts, and no orbit repeats inside a batch.

Those are structural.  The check that actually binds a claim to work is the
double scalar multiplication over `mu = prod_j (1 + s^j)^{n_j}`, and the
authority on it is cairn's own checker, not ours.  Point `--cairn` at a cairn
checkout and this runs `orbit_dp.py verify` against a job document written
from the same generated header; without it that stage is reported as skipped
rather than silently passing.

    python3 test_witness.py
    python3 test_witness.py --cairn ~/aburan28/cairn
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CURVE, INSTANCE, DPW = 41, 2, 12
FAILED = []


def check(label, ok, detail=""):
    print("%-4s %s%s" % ("ok" if ok else "FAIL", label, (" -- " + detail) if detail else ""))
    if not ok:
        FAILED.append(label)
    return ok


def run(cmd, **kw):
    return subprocess.run(cmd, cwd=HERE, capture_output=True, text=True, **kw)


DP_MAGIC_V2 = b"ECC2KDP2"
V1_STRIDE, V2_STRIDE, V2_HEADER = 32, 72, 16


def corpus_framing(path):
    """(base, stride, canon_offset) for whichever format `path` is in."""
    with open(path, "rb") as fh:
        v2 = fh.read(8) == DP_MAGIC_V2
    return (V2_HEADER, V2_STRIDE, 16) if v2 else (0, V1_STRIDE, 8)


def to_v1(src, dst):
    """Rewrite a v2 corpus as the v1 one a counter-less client would have
    written: the same points, with the witness dropped."""
    base, stride, _ = corpus_framing(src)
    raw = open(src, "rb").read()
    if stride == V1_STRIDE:
        open(dst, "wb").write(raw)
        return
    out = bytearray()
    for off in range(base, len(raw), stride):
        rec = raw[off:off + stride]
        out += rec[0:8] + rec[16:40]     # seed, then canon[3]; iters/counts go
    open(dst, "wb").write(bytes(out))


def claim_index(text):
    """{seed: j-vector} over every batch in a witness run's stdout."""
    out = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        for el in json.loads(line)["dps"]:
            out[el["seed"]] = list(el["j"])
    return out


def least_rotation(c, m):
    mask = (1 << m) - 1
    best = c
    for _ in range(m - 1):
        c = ((c << 1) | (c >> (m - 1))) & mask
        best = min(best, c)
    return best


def hex_ok(value):
    return (isinstance(value, str) and value
            and all(ch in "0123456789abcdef" for ch in value)
            and not (len(value) > 1 and value[0] == "0"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cairn", default=os.environ.get("CAIRN_ROOT"),
                    help="a cairn checkout, to cross-verify with its own checker")
    args = ap.parse_args()

    for target in ("cpu", "witness"):
        r = run(["make", "-s", target])
        if not check("build: make %s" % target, r.returncode == 0, r.stderr.strip()[-200:]):
            return 1

    tmp = tempfile.mkdtemp(prefix="witness-")
    corpus = os.path.join(tmp, "dps.bin")
    r = run(["./ecc2k130-cpu", "--curve", str(CURVE), "--instance", str(INSTANCE),
             "--dp-weight", str(DPW), "--dp-file", corpus, "--launches", "400"])
    if not check("client wrote a corpus", os.path.exists(corpus) and os.path.getsize(corpus) > 0,
                 r.stderr.strip()[-200:]):
        return 1
    base, stride, canon_off = corpus_framing(corpus)
    records = (os.path.getsize(corpus) - base) // stride
    check("the corpus frames as %s" % ("v2" if stride == V2_STRIDE else "v1"), records > 0,
          "%d records of %d bytes" % (records, stride))

    job = os.path.join(tmp, "job.json")
    r = run(["python3", "cairn_job.py", "--curve", str(CURVE), "--instance", str(INSTANCE),
             "--dp-weight", str(DPW)], stdout=None)
    check("cairn_job.py writes a job", r.returncode == 0, r.stderr.strip()[-200:])
    open(job, "w").write(r.stdout)
    spec = json.loads(r.stdout)

    # Records early in a corpus are start points that were already
    # distinguished, so their trails are empty and prove nothing about the
    # counts; take the tail, where the walks have actually walked.
    skip = max(0, records - 128)
    r = run(["./build/witness", "--curve", str(CURVE), "--instance", str(INSTANCE),
             "--dp-weight", str(DPW), "--corpus", corpus, "--job", job,
             "--skip", str(skip), "--max", "96", "--batch", "32", "--quiet"])
    if not check("witness ran", r.returncode == 0, r.stderr.strip()[-300:]):
        return 1
    lines = [l for l in r.stdout.splitlines() if l.strip()]

    check("emitted at least one batch", len(lines) >= 1, "%d batches" % len(lines))

    m, maxb = spec["m"], spec["max_batch"]
    total, walked, structural, seen = 0, 0, True, set()
    for line in lines:
        art = json.loads(line)
        if set(art) != {"dps"} or not 1 <= len(art["dps"]) <= maxb:
            structural = False
            continue
        batch_names = set()
        for el in art["dps"]:
            total += 1
            if set(el) != {"x", "seed", "j"}:
                structural = False; continue
            if not hex_ok(el["x"]) or not hex_ok(el["seed"]):
                structural = False; continue
            x = int(el["x"], 16)
            if x >> m or bin(x).count("1") > spec["dp_max_weight"]:
                structural = False; continue
            if x != least_rotation(x, m):
                structural = False; continue
            if not isinstance(el["j"], list) or len(el["j"]) != spec["j_count"]:
                structural = False; continue
            if any(not isinstance(n, int) or isinstance(n, bool) or n < 0 for n in el["j"]):
                structural = False; continue
            if el["x"] in batch_names:
                structural = False; continue
            batch_names.add(el["x"])
            seen.add(el["seed"])
            walked += sum(el["j"])

    check("artifacts are well formed", structural, "%d elements" % total)
    check("every batch is within the job's cap", all(len(json.loads(l)["dps"]) <= maxb for l in lines))
    check("seeds are distinct", len(seen) == total, "%d seeds for %d elements" % (len(seen), total))
    check("the tail of the corpus has walked trails", walked > 0,
          "%d steps across %d trails" % (walked, total))

    # ---- the check this whole path exists for ----------------------------
    # The same trails, read the two ways: straight out of the v2 record, and
    # replayed from a v1 record that dropped the counts.  If the counters the
    # device carried were not the branch counts, these disagree.
    if stride == V2_STRIDE:
        v1corpus = os.path.join(tmp, "dps-v1.bin")
        to_v1(corpus, v1corpus)
        check("the derived v1 corpus has the same record count",
              os.path.getsize(v1corpus) // V1_STRIDE == records,
              "%d vs %d" % (os.path.getsize(v1corpus) // V1_STRIDE, records))
        rv1 = run(["./build/witness", "--curve", str(CURVE), "--instance", str(INSTANCE),
                   "--dp-weight", str(DPW), "--corpus", v1corpus, "--job", job,
                   "--skip", str(skip), "--max", "96", "--batch", "32", "--quiet"])
        if check("witness ran on the replayed v1 corpus", rv1.returncode == 0,
                 rv1.stderr.strip()[-300:]):
            carried, replayed = claim_index(r.stdout), claim_index(rv1.stdout)
            check("both corpus formats cover the same seeds",
                  set(carried) == set(replayed),
                  "%d carried, %d replayed" % (len(carried), len(replayed)))
            same = [k for k in carried if k in replayed and carried[k] == replayed[k]]
            check("carried counts equal replayed counts",
                  len(same) == len(carried) and len(carried) > 0,
                  "%d of %d agree" % (len(same), len(carried)))
    else:
        print("skip carried-vs-replayed -- this client was built without the counters")

    # A corpus record whose orbit the replay does not reach is a disagreement
    # about the walk, and nothing after it is worth claiming.
    bad = os.path.join(tmp, "bad.bin")
    raw = bytearray(open(corpus, "rb").read())
    off = base + skip * stride
    raw[off + canon_off] ^= 0xFF      # flip a bit of the recorded orbit
    open(bad, "wb").write(bytes(raw))
    r = run(["./build/witness", "--curve", str(CURVE), "--instance", str(INSTANCE),
             "--dp-weight", str(DPW), "--corpus", bad, "--job", job,
             "--skip", str(skip), "--max", "1", "--quiet"])
    check("a tampered corpus record is refused", r.returncode != 0,
          "exit %d" % r.returncode)

    # A job naming a normal element that is not a conjugate of this basis has
    # no coordinate permutation, so no name it produced would check.
    alien = json.loads(open(job).read())
    alien["nb_generator"] = "3"
    aliens = os.path.join(tmp, "alien.json")
    open(aliens, "w").write(json.dumps(alien))
    r = run(["./build/witness", "--curve", str(CURVE), "--instance", str(INSTANCE),
             "--dp-weight", str(DPW), "--corpus", corpus, "--job", aliens,
             "--skip", str(skip), "--max", "1", "--quiet"])
    check("a job in a foreign basis is refused", r.returncode != 0, "exit %d" % r.returncode)

    # A job whose constants this binary does not walk is refused before it
    # emits anything at all.
    wrong = json.loads(open(job).read())
    wrong["j_count"] = 7
    wrongs = os.path.join(tmp, "wrong.json")
    open(wrongs, "w").write(json.dumps(wrong))
    r = run(["./build/witness", "--curve", str(CURVE), "--instance", str(INSTANCE),
             "--dp-weight", str(DPW), "--corpus", corpus, "--job", wrongs,
             "--skip", str(skip), "--max", "1", "--quiet"])
    check("a job this binary cannot walk is refused", r.returncode != 0, "exit %d" % r.returncode)

    # ---- the check that binds: cairn's own verifier ----------------------
    orbit_dp = None
    if args.cairn:
        cand = os.path.join(args.cairn, "examples", "certicom-ecdlp", "tools", "orbit_dp.py")
        if os.path.exists(cand):
            orbit_dp = cand
    if orbit_dp is None:
        print("skip cairn cross-verification -- pass --cairn <checkout> or set CAIRN_ROOT")
    else:
        # `verify` reads one artifact per file and the emitter writes one batch
        # per line, so split.  With <= max_batch records the two happen to look
        # the same, which is exactly why this was worth writing out.
        batches = []
        for i, line in enumerate(lines):
            bpath = os.path.join(tmp, "batch-%d.json" % i)
            open(bpath, "w").write(line)
            batches.append(bpath)
        results = [run(["python3", orbit_dp, "verify", "--job", job, b]) for b in batches]
        check("cairn's checker verifies every batch",
              bool(results) and all(x.returncode == 0 for x in results),
              "%d of %d batches" % (sum(x.returncode == 0 for x in results), len(results)))
        # and rejects a witness that does not reach the orbit it names
        art = json.loads(lines[0])
        art["dps"][0]["j"] = list(art["dps"][0]["j"])
        art["dps"][0]["j"][0] += 1
        tampered = os.path.join(tmp, "tampered.json")
        open(tampered, "w").write(json.dumps(art))
        r = run(["python3", orbit_dp, "verify", "--job", job, tampered])
        check("cairn's checker rejects a bumped count", r.returncode != 0,
              (r.stdout + r.stderr).strip()[-160:])

    print()
    if FAILED:
        print("FAILED: %s" % ", ".join(FAILED))
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
