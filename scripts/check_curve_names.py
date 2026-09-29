#!/usr/bin/env python3
"""Enforce ICV1 curve names (docs/curves/ICV1.md, AGENTS.md §11).

    python3 scripts/check_curve_names.py                 # lint Markdown and HTML
    python3 scripts/check_curve_names.py --diff BASE     # also lint what a branch adds
    python3 scripts/check_curve_names.py --fix           # rewrite retired names

The lint fails on

1. a retired curve name (`K_0 / GF(2^41)`, `k0n41`, `bench-24bit`,
   `generated-24bit-10935329`, `random-binary-n27-b845462`,
   `E_{0,2}/GF(4) over GF(2^14)`, and their spellings) in tracked Markdown
   or HTML outside a fenced code block;
2. an `icv1-…` slug in tracked Markdown or HTML that the registry does
   not hold;
3. with `--diff BASE`: a retired name on a line the branch adds to a
   Markdown, HTML, Rust, Python or YAML file, or a file the branch adds
   whose name uses a retired curve stem.

Exempt, because they are evidence and are frozen: files whose SHA-256 a
tracked manifest records, and anything under a `results/`, `runs/`,
`raw/`, `evidence/` or `archives/` directory.  A retired name that is
part of a path (`docs/ic/params/k0n41-subgroup.json`) is a file name,
not a curve name, and is left alone.

`--fix` replaces each retired name with the curve's standard name when
it has one (`ECC2K-130`, `sect163k1`), else its slug.  It is idempotent,
so it can be rerun after merging a branch that added prose.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import curve_id as cid  # noqa: E402

REPO = cid.REPO
TEXT = ("*.md", "*.html")
DIFF_TEXT = ("*.md", "*.html", "*.rs", "*.py", "*.yml", "*.yaml")
FROZEN_DIRS = ("results", "runs", "raw", "evidence", "archives")
# Files whose retired names are the subject, not a use: the spec's table
# of retired forms, and the tools that recognise them.
SELF = {"docs/curves/ICV1.md", "docs/curves/README.md", "scripts/check_curve_names.py",
        "scripts/build_curve_registry.py", "scripts/curve_id.py",
        "src/cryptanalysis/curve_id.rs", "tests/curve_id.rs",
        "examples/curve_id_generated.rs"}

SUB = r"(?:_\{?([01])\}?|<sub>([01])</sub>|&#832([01]);|([₀₁]))"
DEG = (r"(?:GF\(\s*2(?:\^\{?(\d+)\}?|<sup>(\d+)</sup>)\s*\)"
       r"|F_\{?2\^\{?(\d+)\}?\}?|F<sub>2</sub><sup>(\d+)</sup>"
       r"|2(?:\^\{?(\d+)\}?|<sup>(\d+)</sup>))")
KOBLITZ = re.compile(r"K" + SUB + r"\s*/\s*" + DEG)
EDGE_L, EDGE_R = r"(?<![A-Za-z0-9_./-])", r"(?![A-Za-z0-9_.-]*[A-Za-z0-9_])(?![./-][A-Za-z0-9])"
STEM = re.compile(EDGE_L + r"(k[01]n\d{1,3}|bench-\d+bit|generated-\d+bit-\d+"
                  r"|random-binary-n\d+-b[0-9a-f]+)" + EDGE_R)
SUBFIELD = re.compile(r"E_\{(\d+),(\d+)\}/GF\((?:2\^(\d+)|(\d+))\)(?:\s+over\s+GF\(2\^(\d+)\))?")
SLUG = re.compile(r"icv1-(?:f2m|fp)\d+-tm?\d+-[0-9a-f]{8}")
RETIRED_FILE = re.compile(r"(?:^|/)(k[01]n\d+|bench-\d+bit|generated-\d+bit-\d+"
                          r"|random-binary-n\d+-b[0-9a-f]+)[-_.]")


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True,
                          check=True).stdout


def tracked(patterns) -> list[str]:
    return [p for p in git("ls-files", "-z", "--", *patterns).split("\0") if p]


def pinned_hashes() -> dict[str, set[str]]:
    """Every 64-hex-digit string in a tracked manifest, seal, audit or
    checksum file, mapped to the directories of the files that record it."""
    out: dict[str, set[str]] = {}
    for rel in tracked(("*manifest*.json", "*seal*.json", "*audit*.json", "*SHA256SUMS*",
                        "*.sha256")):
        try:
            text = (REPO / rel).read_text(errors="ignore")
        except OSError:
            continue
        for h in set(re.findall(r"[0-9a-f]{64}", text)):
            out.setdefault(h, set()).add(str(Path(rel).parent))
    return out


def exempt(rel: str, text: str, pins: dict[str, set[str]]) -> bool:
    """Frozen evidence: under a results-like directory, or pinned by a
    manifest in its own directory or one above it.  A living document an
    audit elsewhere happens to pin (docs/ic/BOUNDARY_TARGETS.md) is not
    exempt; that audit is re-pinned when the document moves."""
    if rel in SELF:
        return True
    if any(part in FROZEN_DIRS for part in Path(rel).parts[:-1]):
        return True
    here = Path(rel).parent
    for d in pins.get(hashlib.sha256(text.encode()).hexdigest(), ()):
        if d == "." or here == Path(d) or Path(d) in here.parents:
            return True
    return False


def fence_mask(text: str, rel: str) -> list[tuple[int, int]]:
    """Spans of fenced code blocks in Markdown: literal transcripts."""
    if not rel.endswith(".md"):
        return []
    spans, start, pos = [], None, 0
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith("```"):
            if start is None:
                start = pos
            else:
                spans.append((start, pos + len(line)))
                start = None
        pos += len(line)
    return spans


def in_spans(i: int, spans) -> bool:
    return any(a <= i < b for a, b in spans)


class Resolver:
    def __init__(self) -> None:
        self.reg = cid.load_registry()
        self.by_alias: dict[str, dict] = {}
        self.slugs = set()
        for c in self.reg["curves"]:
            self.slugs.add(c["slug"])
            for a in c["aliases"] + c["standard_names"] + [c["slug"]]:
                self.by_alias[cid.normalise_alias(a)] = c

    def koblitz(self, a: int, n: int) -> dict | None:
        return self.by_alias.get(cid.normalise_alias(f"K_{a} / GF(2^{n})"))

    def name(self, c: dict) -> str:
        return c["standard_names"][0] if c["standard_names"] else c["slug"]


def find_retired(text: str, rel: str, res: Resolver):
    """(start, end, found, replacement or None) for each retired name."""
    fences = fence_mask(text, rel)
    hits = []
    for m in KOBLITZ.finditer(text):
        if in_spans(m.start(), fences):
            continue
        a = next(g for g in m.group(1, 2, 3, 4) if g)
        a = int("₀₁".index(a)) if a in "₀₁" else int(a)
        n = int(next(g for g in m.group(5, 6, 7, 8, 9, 10) if g))
        c = res.koblitz(a, n)
        hits.append((m.start(), m.end(), m.group(0), res.name(c) if c else None))
    for m in STEM.finditer(text):
        if in_spans(m.start(), fences):
            continue
        tok = m.group(1)
        if tok.startswith("k"):
            c = res.koblitz(int(tok[1]), int(tok[3:]))
        else:
            c = res.by_alias.get(cid.normalise_alias(tok))
        hits.append((m.start(), m.end(), tok, res.name(c) if c else None))
    for m in SUBFIELD.finditer(text):
        if in_spans(m.start(), fences):
            continue
        a, b = m.group(1), m.group(2)
        k = int(m.group(3)) if m.group(3) else {2: 1, 4: 2, 8: 3, 16: 4}.get(int(m.group(4)))
        c = None
        if m.group(5) and k:
            c = res.by_alias.get(cid.normalise_alias(f"E_{{{a},{b}}}/GF(2^{k}) over GF(2^{m.group(5)})"))
        hits.append((m.start(), m.end(), m.group(0), res.name(c) if c else None))
    return sorted(hits)


def lint(res: Resolver, pins: dict[str, set[str]], fix: bool) -> list[str]:
    problems, changed = [], 0
    for rel in tracked(TEXT):
        path = REPO / rel
        text = path.read_text(errors="ignore")
        if exempt(rel, text, pins):
            continue
        hits = find_retired(text, rel, res)
        if fix and hits:
            out, last = [], 0
            for a, b, found, repl in hits:
                if a < last:
                    continue
                if repl is None:
                    problems.append(f"{rel}: cannot rewrite {found!r}: not in the registry")
                    continue
                out.append(text[last:a])
                out.append(repl)
                last = b
            out.append(text[last:])
            new = "".join(out)
            if new != text:
                path.write_text(new)
                changed += 1
            text = new
            hits = find_retired(text, rel, res)
        for a, _, found, repl in hits:
            line = text.count("\n", 0, a) + 1
            hint = f" (write {repl})" if repl else " (not in the registry)"
            problems.append(f"{rel}:{line}: retired curve name {found!r}{hint}")
        for m in SLUG.finditer(text):
            if m.group(0) not in res.slugs:
                line = text.count("\n", 0, m.start()) + 1
                problems.append(f"{rel}:{line}: {m.group(0)} is not in docs/curves/registry.json")
    if fix:
        print(f"rewrote {changed} files", file=sys.stderr)
    return problems


def lint_diff(base: str, res: Resolver) -> list[str]:
    problems = []
    for rel in git("diff", "--name-only", "--diff-filter=A", f"{base}...HEAD").split():
        if RETIRED_FILE.search(rel):
            problems.append(f"{rel}: new file named by a retired curve stem; name it by the slug")
    diff = git("diff", "-U0", f"{base}...HEAD", "--", *DIFF_TEXT)
    rel, line = None, 0
    for raw in diff.splitlines():
        if raw.startswith("+++ "):
            rel = raw[6:] if raw.startswith("+++ b/") else None
        elif raw.startswith("@@"):
            line = int(re.search(r"\+(\d+)", raw).group(1))
        elif raw.startswith("+") and rel and rel not in SELF:
            for _, _, found, repl in find_retired(raw[1:], "added.txt", res):
                hint = f" (write {repl})" if repl else ""
                problems.append(f"{rel}:{line}: added line uses retired curve name {found!r}{hint}")
            line += 1
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fix", action="store_true", help="rewrite retired names in place")
    ap.add_argument("--diff", metavar="BASE", help="also lint what HEAD adds over BASE")
    args = ap.parse_args()
    res = Resolver()
    problems = lint(res, pinned_hashes(), args.fix)
    if args.diff:
        problems += lint_diff(args.diff, res)
    for p in problems:
        print(p)
    print(f"curve names: {len(problems)} problem(s)", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
