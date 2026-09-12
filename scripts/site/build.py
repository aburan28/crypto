#!/usr/bin/env python3
"""Assemble the published GitHub Pages site for aburan28/crypto.

The hourly ECC2K-130 workflow uploads the directory this script writes.
Everything it publishes already exists in the repository; the script only
copies and lays out, so the repository file stays canonical (AGENTS.md §7)
and no published page is ever the only copy of a figure.

Layout, and why each path is where it is:

    /                     landing page
    /assets/site.css      landing-page styles
    /scoreboard/          docs/index-calculus-scoreboard.html, the cost ledger
    /status/              the ECC2K-130 distinguished-point dashboard
    /status/status.json   snapshot, next to the page that reads it
    /status/history.json  hourly history, likewise
    /status.json          same bytes at the root
    /history.json         same bytes at the root

The two root-level JSON copies are a compatibility contract, not a
convenience: the publish workflow reads the previous
https://aburan28.github.io/crypto/history.json to merge history forward, so
moving that file loses every published snapshot before the move. The copies
under /status/ let the dashboard use relative fetches and so keep working
when opened straight from the working tree.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

BASE_URL = "https://aburan28.github.io/crypto"

# (source relative to repo root, destination relative to the site root).
PAGES = (
    ("docs/site/index.html", "index.html"),
    ("docs/site/404.html", "404.html"),
    ("docs/site/assets/site.css", "assets/site.css"),
    ("docs/index-calculus-scoreboard.html", "scoreboard/index.html"),
    ("docs/ecc2k130-status/index.html", "status/index.html"),
    ("docs/ecc2k130-status/style.css", "status/style.css"),
)

# Snapshot data, published twice: beside the dashboard and at the root.
DATA = (
    ("docs/ecc2k130-status/status.json", ("status/status.json", "status.json")),
    ("docs/ecc2k130-status/history.json", ("status/history.json", "history.json")),
)

# Pages worth listing for crawlers. Data files and the 404 stay out.
SITEMAP = ("/", "/scoreboard/", "/status/")


def copy(src_rel: str, dest_rel: str, out_dir: str, root: str = ROOT) -> str:
    src = os.path.join(root, src_rel)
    if not os.path.exists(src):
        raise SystemExit(f"missing site source: {src_rel}")
    dest = os.path.join(out_dir, dest_rel)
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    shutil.copyfile(src, dest)
    return dest


def render_robots() -> str:
    return "\n".join(
        (
            "User-agent: *",
            "Allow: /",
            f"Sitemap: {BASE_URL}/sitemap.xml",
            "",
        )
    )


def render_sitemap(lastmod: str) -> str:
    lines = ['<?xml version="1.0" encoding="UTF-8"?>', '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
    for path in SITEMAP:
        lines += ["  <url>", f"    <loc>{BASE_URL}{path}</loc>", f"    <lastmod>{lastmod}</lastmod>", "  </url>"]
    lines += ["</urlset>", ""]
    return "\n".join(lines)


def write(text: str, dest_rel: str, out_dir: str) -> str:
    dest = os.path.join(out_dir, dest_rel)
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    with open(dest, "w", encoding="utf-8") as fh:
        fh.write(text)
    return dest


def build(out_dir: str, root: str = ROOT, lastmod: str | None = None) -> list[str]:
    """Write the site into out_dir and return the paths published, sorted."""
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    written = [copy(src, dest, out_dir, root) for src, dest in PAGES]
    for src, dests in DATA:
        written += [copy(src, dest, out_dir, root) for dest in dests]

    stamp = lastmod or dt.datetime.now(dt.timezone.utc).date().isoformat()
    written.append(write(render_robots(), "robots.txt", out_dir))
    written.append(write(render_sitemap(stamp), "sitemap.xml", out_dir))
    # Pages is served by Actions here, so Jekyll never runs; the marker keeps
    # that true if the source is ever switched back to a branch.
    written.append(write("", ".nojekyll", out_dir))
    return sorted(os.path.relpath(path, out_dir) for path in written)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Assemble the GitHub Pages site.")
    parser.add_argument("--out", default="_site", help="output directory (recreated)")
    parser.add_argument("--root", default=ROOT, help="repository root to read sources from")
    parser.add_argument("--lastmod", default=None, help="sitemap lastmod date (default: today, UTC)")
    args = parser.parse_args(argv)

    published = build(os.path.abspath(args.out), os.path.abspath(args.root), args.lastmod)
    for path in published:
        print(path)
    print(f"{len(published)} files -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
