#!/usr/bin/env python3
"""Replay a cell with the tracked binary and compare its counted units with
the frozen raw results: label, found, calls, tame, wild, XOR words and the
tame depths must all agree, target by target.

    python3 replay_check.py --n 7 --encodings SUB,FC --targets 4
"""
import argparse, json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "hamming_pdp", "target", "release", "hamming_pdp")
KEYS = ["label", "found", "agree", "calls", "tame", "wild", "budget_calls", "matrix_cap_hits", "xor_words", "tame_depths", "vars", "gens", "exhaustive_ordered_pairs"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--encodings", default="SUB,FC")
    ap.add_argument("--targets", type=int, default=4)
    ap.add_argument("--tag", default="main")
    a = ap.parse_args()
    bad = 0
    for enc in a.encodings.split(","):
        frozen_path = os.path.join(HERE, "results", a.tag, f"n{a.n}_{enc}.jsonl")
        frozen = {}
        manifest = None
        with open(frozen_path) as fh:
            for line in fh:
                r = json.loads(line)
                if r["kind"] == "manifest":
                    manifest = r
                else:
                    frozen[r["seed"]] = r
        with tempfile.NamedTemporaryFile(suffix=".jsonl", delete=False) as tf:
            out = tf.name
        cmd = [BIN, "--n", str(a.n), "--w", str(manifest["w"]), "--l", str(manifest["l"]), "--targets", str(a.targets),
               "--encodings", enc, "--budget-log2", str(manifest["budget_xor"].bit_length() - 1),
               "--matrix-cap-log2", str(manifest["matrix_cap_words"].bit_length() - 1), "--max-calls", str(manifest["max_calls"]), "--out", out]
        subprocess.check_call(cmd, stderr=subprocess.DEVNULL)
        with open(out) as fh:
            for line in fh:
                r = json.loads(line)
                if r["kind"] == "manifest":
                    for k in ["irr", "normal_alpha", "order", "wt_x_digest", "sub_x_digest", "sub_basis"]:
                        if r[k] != manifest[k]:
                            print(f"MANIFEST {enc} {k}: {r[k]} != {manifest[k]}"); bad += 1
                    continue
                f = frozen[r["seed"]]
                for k in KEYS:
                    if r[k] != f[k]:
                        print(f"DIFF n={a.n} {enc} seed={r['seed']} {k}: replay {r[k]} frozen {f[k]}"); bad += 1
                print(f"n={a.n} {enc} seed={r['seed']} calls={r['calls']} xor={r['xor_words']} {'ok' if all(r[k] == f[k] for k in KEYS) else 'DIFF'}")
    print("replay differences:", bad)
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
