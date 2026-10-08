#!/usr/bin/env python3
"""Write the panel index on the current-state page from the scoreboard.

Every top-level panel of docs/index-calculus-scoreboard.html carries a
one-line summary under its title.  This script lists those panels, in ledger
order, between the ``<!-- BEGIN panel-index -->`` and ``<!-- END panel-index -->``
markers of docs/ic-current-state.html: title linked to the panel, its class
chip, and its summary line.  The index cites the ledger and computes nothing.

    python3 scripts/site/panel_index.py          # rewrite the index
    python3 scripts/site/panel_index.py --check  # fail if it is stale
"""

from __future__ import annotations

import argparse
import html
import os
import re
import sys
from html.parser import HTMLParser

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SCOREBOARD = os.path.join(ROOT, "docs", "index-calculus-scoreboard.html")
PAGE = os.path.join(ROOT, "docs", "ic-current-state.html")
BEGIN = "<!-- BEGIN panel-index -->"
END = "<!-- END panel-index -->"

VOID = {"area", "base", "br", "col", "embed", "hr", "img", "input",
        "link", "meta", "param", "source", "track", "wbr"}


class Walker(HTMLParser):
    """Collect each top-level panel's id, title, class chip and summary."""

    def __init__(self):
        super().__init__()
        self.stack = []        # (tag, classes)
        self.panels = []
        self.cur = None        # the panel being read
        self.target = None     # "title", "chip" or "summary" while reading text

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        classes = (attrs.get("class") or "").split()
        panel = "panel" in classes or (
            tag == "section" and bool(attrs.get("id") or attrs.get("aria-labelledby")))
        outside = not any(t in ("header", "aside") or "notes" in c for t, c in self.stack)
        if panel and self.cur is None and outside:
            self.cur = {"id": attrs.get("id") or "", "depth": len(self.stack),
                        "title": [], "chips": [], "summary": [], "heading_id": None}
        if self.cur is not None:
            if tag in ("h2", "h3") and self.target is None and not self.cur["title"]:
                self.target = "title"
                self.cur["heading_id"] = attrs.get("id")
            elif self.target == "title" and "chip" in classes:
                self.cur["chips"].append([])
                self.target = "chip"
            elif tag == "p" and "panel-summary" in classes:
                self.target = "summary"
        if tag not in VOID:
            self.stack.append((tag, classes))

    def handle_data(self, data):
        if self.cur is None or self.target is None:
            return
        if self.target == "title":
            self.cur["title"].append(data)
        elif self.target == "chip":
            self.cur["chips"][-1].append(data)
        elif self.target == "summary":
            self.cur["summary"].append(data)

    def handle_endtag(self, tag):
        if tag in VOID:
            return
        if self.cur is not None:
            if self.target == "chip" and tag == "span":
                self.target = "title"
            elif self.target == "title" and tag in ("h2", "h3"):
                self.target = "closed"
            elif self.target == "summary" and tag == "p":
                self.target = "closed"
        for i in range(len(self.stack) - 1, -1, -1):
            if self.stack[i][0] == tag:
                del self.stack[i:]
                break
        if self.cur is not None and len(self.stack) == self.cur["depth"]:
            self.panels.append(self.finish(self.cur))
            self.cur = None
            self.target = None

    @staticmethod
    def finish(cur):
        squash = lambda parts: re.sub(r"\s+", " ", "".join(parts)).strip()
        return {
            "id": cur["id"] or cur["heading_id"] or "",
            "title": squash(cur["title"]),
            "chips": [squash(c) for c in cur["chips"]],
            "summary": squash(cur["summary"]),
        }


def panels(markup):
    walker = Walker()
    walker.feed(markup)
    walker.close()
    return walker.panels


def render(items):
    out = [BEGIN,
           "  <details>",
           "    <summary>Show all %d panels, one line each</summary>" % len(items),
           '    <ol class="index">']
    for p in items:
        if not p["id"]:
            raise ValueError("a scoreboard panel has no id to link to: %r" % p["title"])
        if not p["summary"]:
            raise ValueError("scoreboard panel %s has no panel-summary line" % p["id"])
        chips = "".join(' <span class="chip">%s</span>' % html.escape(c) for c in p["chips"])
        out.append('      <li><a href="./#%s">%s</a>%s<br>%s</li>' % (
            html.escape(p["id"], quote=True), html.escape(p["title"]), chips,
            html.escape(p["summary"])))
    out += ["    </ol>", "  </details>", END]
    return "\n".join(out)


def splice(page, section):
    if page.count(BEGIN) != 1 or page.count(END) != 1:
        raise ValueError("docs/ic-current-state.html must carry one panel-index region")
    start = page.index(BEGIN)
    end = page.index(END) + len(END)
    return page[:start] + section + page[end:]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="verify instead of writing")
    args = parser.parse_args(argv)
    with open(SCOREBOARD, encoding="utf-8") as fh:
        items = panels(fh.read())
    with open(PAGE, encoding="utf-8") as fh:
        page = fh.read()
    new = splice(page, render(items))
    if args.check:
        if new != page:
            print("panel index on docs/ic-current-state.html is stale; "
                  "run scripts/site/panel_index.py", file=sys.stderr)
            return 1
        print("panel index: %d panels, current" % len(items))
        return 0
    with open(PAGE, "w", encoding="utf-8") as fh:
        fh.write(new)
    print("panel index: %d panels written" % len(items))
    return 0


if __name__ == "__main__":
    sys.exit(main())
