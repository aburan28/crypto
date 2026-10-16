#!/usr/bin/env python3
"""Run msolve F4 on exported chained-system draws and record D_F4.

usage: python3 -I run_f4.py OUT.jsonl JOB [JOB ...]
JOB = m:n:ell:seed:k[:sat|:nocount]   (k = unsat index; with ':sat', the k-th
draw regardless of satisfiability, for the instrument-correctness arm E; with
':nocount', the k-th draw with no exact count, msolve's GB deciding)
     | rand:FILE.ms:label     (a pre-built msolve file, arm D)
Each job: export (ic_next10_06_chain_export), msolve -t 1 -g 2 -v 2 under
prlimit --as=3e9 and a 2700 s timeout, parse the per-round degrees.
"""
import json, os, re, resource, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../.."))
EXPORT = os.path.join(ROOT, "target/release/examples/ic_next10_06_chain_export")
WORK = os.path.join(HERE, "runs", "work")
CAP = "3000000000"
# MSOLVE_EXTRA, e.g. "-u 1": rebuild the hash table every step (memory only;
# the computation and its degrees are unchanged).
EXTRA = os.environ.get("MSOLVE_EXTRA", "").split()
TIMEOUT = int(os.environ.get("F4_TIMEOUT", "2700"))
ROUND = re.compile(r"^\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+) x (\d+)\s")


def run_msolve(ms, tag):
    log = os.path.join(WORK, tag + ".log")
    out = os.path.join(WORK, tag + ".out")
    env = dict(os.environ, OMP_NUM_THREADS="1", RAYON_NUM_THREADS="1")
    t0 = time.time()
    with open(log, "w") as lf:
        p = subprocess.Popen(["prlimit", "--as=" + CAP, "--", "nice", "msolve", "-t", "1",
                              "-g", "2", "-v", "2"] + EXTRA + ["-f", ms, "-o", out],
                             stdout=lf, stderr=subprocess.STDOUT, env=env)
        while True:
            pid, status, ru = os.wait4(p.pid, os.WNOHANG)
            if pid:
                code = os.waitstatus_to_exitcode(status)
                rss = ru.ru_maxrss  # this child's peak RSS, kB
                cpu = ru.ru_utime + ru.ru_stime
                break
            if time.time() - t0 > TIMEOUT:
                p.kill()
                _, status, ru = os.wait4(p.pid, 0)
                code, rss = "timeout", ru.ru_maxrss
                cpu = ru.ru_utime + ru.ru_stime
                break
            time.sleep(0.05)
    secs = time.time() - t0
    degs, maxmat = [], (0, 0)
    for line in open(log, errors="replace"):
        mm = ROUND.match(line)
        if mm:
            degs.append(int(mm.group(1)))
            r, c = int(mm.group(4)), int(mm.group(5))
            if r * c > maxmat[0] * maxmat[1]:
                maxmat = (r, c)
    gb_one = None
    if os.path.exists(out):
        txt = open(out).read()
        body = txt.split("#---")[-1].strip()
        gb_one = body.startswith("[1]")
        size = re.search(r"#length of basis:\s+(\d+)", txt)
        gb_len = int(size.group(1)) if size else None
    else:
        gb_len = None
    return {"exit": code, "secs": round(secs, 3), "cpu_secs": round(cpu, 3), "peak_rss_kb": rss,
            "rounds": len(degs), "degrees": degs, "d_f4": max(degs) if degs else None,
            "max_matrix": list(maxmat), "gb_is_one": gb_one, "gb_len": gb_len}


def main():
    os.makedirs(WORK, exist_ok=True)
    outpath = sys.argv[1]
    for job in sys.argv[2:]:
        parts = job.split(":")
        if parts[0] == "rand":
            rec = {"job": job, "label": parts[2]}
            ms, tag = parts[1], "rand-" + parts[2]
        else:
            m, n, ell, seed, k = parts[:5]
            sat = len(parts) > 5 and parts[5] == "sat"
            tag = f"m{m}-n{n}-l{ell}-s{seed}-{'d' if sat else 'u'}{k}"
            ms = os.path.join(WORK, tag + ".ms")
            cmd = [EXPORT, "--m", m, "--n", n, "--ell", ell, "--seed", seed,
                   "--unsat-index", k, "--out", ms, "--max-draws", "20000"]
            if sat:
                cmd.append("--sat-ok")
            if len(parts) > 5 and parts[5] == "nocount":
                cmd += ["--sat-ok", "--no-count"]
                tag = tag.replace("-u", "-nc")
            t0 = time.time()
            p = subprocess.run(cmd, capture_output=True, text=True)
            if p.returncode != 0:
                print(json.dumps({"job": job, "export_failed": p.stderr.strip()}), flush=True)
                with open(outpath, "a") as f:
                    f.write(json.dumps({"job": job, "export_failed": p.stderr.strip()}) + "\n")
                continue
            rec = {"job": job, "export": json.loads(p.stdout.strip()),
                   "export_secs": round(time.time() - t0, 3)}
        rec.update(run_msolve(ms, tag))
        rec["msolve_extra"] = " ".join(EXTRA)
        if parts[0] != "rand":
            os.remove(ms)  # reproducible from the exporter; keeps raw data small
        line = json.dumps(rec)
        print(line, flush=True)
        with open(outpath, "a") as f:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
