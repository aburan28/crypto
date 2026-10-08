#!/usr/bin/env python3
"""Refresh the Markdown tables between `<!-- TABLE:name -->` and
`<!-- /TABLE -->` markers in RESULT.md and the research note from the
frozen summary (render_tables.py) and the budget-sweep table.

    python3 fill_tables.py
"""
import os, re, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FILES = [
    os.path.join(HERE, "RESULT.md"),
    os.path.join(REPO, "research", "notes", "index-calculus", "RESEARCH_HAMMING_IDEAL_PDP.md"),
]

def table(name):
    if name == "tau":
        return open(os.path.join(HERE, "results", "tau", "TABLE.md")).read().strip()
    return subprocess.check_output(["python3", os.path.join(HERE, "render_tables.py"), name], text=True).strip()

def main():
    for path in FILES:
        s = open(path).read()
        def repl(m):
            return f"<!-- TABLE:{m.group(1)} -->\n{table(m.group(1))}\n<!-- /TABLE -->"
        s2 = re.sub(r"<!-- TABLE:(\w+) -->\n.*?<!-- /TABLE -->", repl, s, flags=re.S)
        open(path, "w").write(s2)
        print("filled", os.path.relpath(path, REPO))

if __name__ == "__main__":
    main()
