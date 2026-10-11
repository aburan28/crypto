#!/usr/bin/env python3
"""Fail when a workflow's path filter misses a file the leaderboard builder reads.

    python3 scripts/check_ic_leaderboard_paths.py [workflow.yml ...]
                                    # default: .github/workflows/ic-leaderboard.yml

`scripts/build_ic_leaderboard.py --check` is only as good as the pull requests
that run it.  A workflow triggered by `paths:` runs on a change to a file the
list names and on nothing else, so a frozen file the builder reads but the
list omits can change the page's numbers with no check run, and the page goes
stale until some unrelated change happens to trigger one.  The list is kept by
hand and the builder's inputs grow with every round, so this runs the builder
under a read recorder and tests every file it opened against the workflow's
`pull_request` and `push` lists.  The builder is itself on both lists, so the
change that adds an input also runs this check.

Stdlib only; the workflow's filter is read as text.  Only the glob forms the
lists use are modelled: literal text, `*` (stops at `/`) and `**` (does not).
GitHub also gives `?`, `+`, `[...]` and a leading `!` meanings of their own;
a pattern using one is refused rather than matched wrongly.
"""
from __future__ import annotations

import builtins
import collections
import io
import os
import re
import runpy
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
BUILDER = "scripts/build_ic_leaderboard.py"
DEFAULT = [".github/workflows/ic-leaderboard.yml"]
EVENTS = ("pull_request", "push")
# What the builder writes.  Editing them by hand is caught by --check on any run;
# whether a change to one of them triggers a run is the workflow's own business.
OUTPUTS = {"docs/ic/leaderboard.json", "docs/ic/LEADERBOARD.md", "docs/ic-leaderboard.html"}


def builder_reads() -> list[str]:
    """Every file under the repository that `build_ic_leaderboard.py --check` opens."""
    opened: set[str] = set()
    real_open = builtins.open
    real_bytes, real_text = Path.read_bytes, Path.read_text

    def note(p) -> None:
        if isinstance(p, (str, bytes, os.PathLike)):
            opened.add(os.path.realpath(p))

    def spy_open(f, *a, **k):
        note(f)
        return real_open(f, *a, **k)

    def spy_bytes(self, *a, **k):
        note(self)
        return real_bytes(self, *a, **k)

    def spy_text(self, *a, **k):
        note(self)
        return real_text(self, *a, **k)

    builtins.open, Path.read_bytes, Path.read_text = spy_open, spy_bytes, spy_text
    argv, cwd, path = sys.argv, os.getcwd(), list(sys.path)
    try:
        os.chdir(REPO)
        sys.path.insert(0, str(REPO / "scripts"))
        sys.argv = [BUILDER, "--check"]
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            try:
                runpy.run_path(str(REPO / BUILDER), run_name="__main__")
            except SystemExit:
                pass  # a stale output is --check's business, not this script's
    finally:
        builtins.open, Path.read_bytes, Path.read_text = real_open, real_bytes, real_text
        sys.argv, sys.path[:] = argv, path
        os.chdir(cwd)
    root = str(REPO) + os.sep
    rel = sorted(f[len(root):] for f in opened if f.startswith(root))
    return [r for r in rel if not r.startswith(".git/") and r not in OUTPUTS]


def glob_regex(pat: str) -> re.Pattern:
    if re.search(r"[?+\[\]!]", pat):
        raise SystemExit(f"path pattern {pat!r} uses a GitHub glob form this check does not model "
                         "(? + [] !); list the files with `*` and `**` instead")
    out, i = "", 0
    while i < len(pat):
        if pat.startswith("**/", i):
            out, i = out + "(?:.*/)?", i + 3
        elif pat.startswith("**", i):
            out, i = out + ".*", i + 2
        elif pat[i] == "*":
            out, i = out + "[^/]*", i + 1
        else:
            out, i = out + re.escape(pat[i]), i + 1
    return re.compile("^" + out + "$")


def event_paths(workflow: Path, event: str) -> list[str]:
    """The `paths:` list under `on: <event>:` of a workflow, read as text."""
    patterns: list[str] = []
    in_event = in_paths = False
    for line in workflow.read_text().splitlines():
        if re.match(r"^  \w[\w-]*:", line):  # a key of `on:`
            in_event, in_paths = line.strip().startswith(event + ":"), False
            continue
        if not in_event:
            continue
        if re.match(r"^    paths-ignore:", line):
            raise SystemExit(f"{workflow.name}: {event} uses paths-ignore; this check does not model it")
        if re.match(r"^    paths:", line):
            in_paths = True
        elif re.match(r"^    \S", line):
            in_paths = False
        elif in_paths:
            m = re.match(r"^\s+-\s+(?:'([^']+)'|\"([^\"]+)\"|([^\s#]+))", line)
            if m:
                patterns.append(next(g for g in m.groups() if g))
    return patterns


def label(path: str) -> str:
    """Group the many files of one run directory on one line."""
    return re.sub(r"(/runs[^/]*(?:/main)?)/.*", r"\1/…", path)


def main() -> int:
    workflows = [Path(a) if Path(a).is_absolute() else REPO / a for a in (sys.argv[1:] or DEFAULT)]
    reads = builder_reads()
    if not reads:
        print("the builder read no files; the recorder is broken", file=sys.stderr)
        return 2
    bad = 0
    for wf in workflows:
        rel = wf.relative_to(REPO) if wf.is_relative_to(REPO) else wf
        for event in EVENTS:
            patterns = event_paths(wf, event)
            if not patterns:
                print(f"{rel}: no `paths:` list for {event}", file=sys.stderr)
                bad += 1
                continue
            rx = [glob_regex(p) for p in patterns]
            missed = [r for r in reads if not any(p.match(r) for p in rx)]
            if BUILDER not in patterns:
                missed.append(BUILDER)
            if missed:
                bad += 1
                groups = collections.Counter(label(m) for m in missed)
                print(f"{rel} [{event}] does not cover {len(missed)} of {len(reads)} files the builder "
                      f"reads; add to its `paths:`:", file=sys.stderr)
                for g, n in sorted(groups.items()):
                    print(f"  {g}" + (f"   ({n} files)" if n > 1 else ""), file=sys.stderr)
            else:
                print(f"{rel} [{event}] covers all {len(reads)} files the builder reads")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
