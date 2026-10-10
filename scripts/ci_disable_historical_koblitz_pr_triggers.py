#!/usr/bin/env python3
"""Rewrite koblitz-stage*.yml triggers to workflow_dispatch only.

Historical stage workflows shared path filters such as
``examples/koblitz_rank_fixture.rs``. One PR touching that file used to
enqueue dozens of stage checks and starve GitHub-hosted runners across
every repository on the account. Re-run this script after adding a new
``koblitz-stage*.yml`` if it ships with a ``pull_request`` trigger.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
WORKFLOWS = ROOT / ".github" / "workflows"
NOTE = (
    "  # Historical stage control. Manual only: shared fixture path\n"
    "  # filters used to enqueue dozens of these on one PR and starve\n"
    "  # GitHub-hosted runners account-wide.\n"
)


def rewrite(text: str) -> str:
    match = re.search(r"(?m)^on:\n", text)
    if not match:
        raise ValueError("no on: block")
    start = match.end()
    lines = text[start:].splitlines(keepends=True)
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() == "":
            j = i + 1
            while j < len(lines) and lines[j].strip() == "":
                j += 1
            if j < len(lines) and not lines[j].startswith((" ", "\t")) and lines[j].strip():
                break
            i += 1
            continue
        if not line.startswith(("  ", "\t")):
            break
        i += 1
    return text[: match.start()] + "on:\n" + NOTE + "  workflow_dispatch:\n" + "".join(lines[i:])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="exit 1 if any koblitz-stage workflow still has a non-dispatch trigger",
    )
    args = parser.parse_args(argv)
    paths = sorted(WORKFLOWS.glob("koblitz-stage*.yml"))
    if not paths:
        print("no koblitz-stage workflows at", WORKFLOWS, file=sys.stderr)
        return 1
    changed = 0
    bad = []
    for path in paths:
        original = path.read_text(encoding="utf-8")
        updated = rewrite(original)
        if args.check:
            if updated != original:
                bad.append(path.name)
            continue
        if updated != original:
            path.write_text(updated, encoding="utf-8")
            changed += 1
            print("rewrote", path.relative_to(ROOT))
    if args.check:
        if bad:
            print("still have PR (or other) triggers:", ", ".join(bad), file=sys.stderr)
            return 1
        print("%d koblitz-stage workflows are workflow_dispatch-only" % len(paths))
        return 0
    print("rewrote %d of %d" % (changed, len(paths)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
