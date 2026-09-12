#!/usr/bin/env python3
"""Offline tests for the GitHub Pages site assembly."""

from __future__ import annotations

import os
import re
import tempfile
import unittest

from build import BASE_URL, build

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# GitHub Pages serves this project under /crypto/, so 404.html (which is
# served for arbitrary depths and cannot use relative links) hardcodes it.
PROJECT_PREFIX = "/crypto/"


def read(path, mode="r"):
    kwargs = {"encoding": "utf-8"} if mode == "r" else {}
    with open(path, mode, **kwargs) as fh:
        return fh.read()


class BuildTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.out = os.path.join(cls._tmp.name, "_site")
        cls.published = build(cls.out, ROOT, lastmod="2026-01-01")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_publishes_every_expected_path(self):
        for path in (
            "index.html",
            "404.html",
            "assets/site.css",
            "scoreboard/index.html",
            "status/index.html",
            "status/style.css",
            "status/status.json",
            "status/history.json",
            "status.json",
            "history.json",
            "robots.txt",
            "sitemap.xml",
            ".nojekyll",
        ):
            self.assertIn(path, self.published, path)
            self.assertTrue(os.path.exists(os.path.join(self.out, path)), path)

    def test_history_json_stays_at_the_site_root(self):
        # The publish job reads the previous history from
        # <BASE_URL>/history.json to merge snapshots forward. Moving that file
        # silently discards every published snapshot, so pin the contract here
        # and in the workflow together.
        workflow = read(os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml"))
        self.assertIn("PAGES_HISTORY_URL: %s/history.json" % BASE_URL, workflow)
        self.assertTrue(os.path.exists(os.path.join(self.out, "history.json")))

    def test_snapshot_data_is_byte_identical_in_both_locations(self):
        for name in ("status.json", "history.json"):
            source = read(os.path.join(ROOT, "docs", "ecc2k130-status", name), "rb")
            self.assertEqual(read(os.path.join(self.out, name), "rb"), source, name)
            self.assertEqual(read(os.path.join(self.out, "status", name), "rb"), source, name)

    def test_scoreboard_is_a_copy_of_the_canonical_repository_file(self):
        # AGENTS.md §7: the repository file is canonical and a published copy
        # is only ever a republish of it.
        source = read(os.path.join(ROOT, "docs", "index-calculus-scoreboard.html"), "rb")
        self.assertEqual(read(os.path.join(self.out, "scoreboard", "index.html"), "rb"), source)

    def test_dashboard_fetches_resolve_next_to_the_published_page(self):
        page = read(os.path.join(self.out, "status", "index.html"))
        for target in re.findall(r'fetch\("\./([^"]+)"', page):
            self.assertTrue(
                os.path.exists(os.path.join(self.out, "status", target)),
                "dashboard fetches %s, which the build does not publish beside it" % target,
            )

    def test_pages_do_not_render_fetched_data_through_innerhtml(self):
        # The dashboard builds rows with textContent because worker_id is
        # self-reported by the walkers; keep it that way.
        for rel in ("index.html", "status/index.html"):
            page = read(os.path.join(self.out, rel))
            self.assertNotIn("innerHTML", page, rel)

    def test_internal_links_resolve(self):
        missing = []
        for rel in ("index.html", "404.html", "status/index.html", "scoreboard/index.html"):
            page = read(os.path.join(self.out, rel))
            base = os.path.dirname(rel)
            for href in re.findall(r'(?:href|src)="([^"]+)"', page):
                if href.startswith(("http://", "https://", "mailto:", "data:", "#", "//")):
                    continue
                target = href.split("#", 1)[0].split("?", 1)[0]
                if not target:
                    continue
                if target.startswith(PROJECT_PREFIX):
                    resolved = target[len(PROJECT_PREFIX):]
                elif target.startswith("/"):
                    missing.append("%s -> %s (root-relative link ignores the /crypto/ prefix)" % (rel, href))
                    continue
                else:
                    resolved = os.path.normpath(os.path.join(base, target))
                candidate = os.path.join(self.out, resolved) if resolved not in ("", ".") else self.out
                if os.path.isdir(candidate):
                    candidate = os.path.join(candidate, "index.html")
                if not os.path.exists(candidate):
                    missing.append("%s -> %s" % (rel, href))
        self.assertEqual(missing, [], "unresolvable internal links: %s" % missing)

    def test_every_page_declares_title_viewport_and_description(self):
        for rel in ("index.html", "404.html", "status/index.html", "scoreboard/index.html"):
            page = read(os.path.join(self.out, rel))
            self.assertRegex(page, r"<title>[^<]+</title>", rel)
            self.assertIn('name="viewport"', page, rel)
            self.assertIn('name="description"', page, rel)

    def test_sitemap_and_robots_point_at_the_published_urls(self):
        sitemap = read(os.path.join(self.out, "sitemap.xml"))
        for path in ("/", "/scoreboard/", "/status/"):
            self.assertIn("<loc>%s%s</loc>" % (BASE_URL, path), sitemap)
        self.assertIn("<lastmod>2026-01-01</lastmod>", sitemap)
        self.assertIn("Sitemap: %s/sitemap.xml" % BASE_URL, read(os.path.join(self.out, "robots.txt")))

    def test_build_is_idempotent(self):
        again = build(self.out, ROOT, lastmod="2026-01-01")
        self.assertEqual(again, self.published)

    def test_missing_source_fails_loudly(self):
        with tempfile.TemporaryDirectory() as empty:
            with self.assertRaises(SystemExit):
                build(os.path.join(empty, "_site"), empty)


class WorkflowTests(unittest.TestCase):
    def test_publish_job_uploads_the_assembled_site(self):
        workflow = read(os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml"))
        self.assertIn("python3 scripts/site/build.py --out _site", workflow)
        self.assertIn("path: _site", workflow)
        self.assertNotIn("path: docs/ecc2k130-status", workflow)

    def test_site_sources_trigger_the_workflow(self):
        workflow = read(os.path.join(ROOT, ".github", "workflows", "ecc2k130-status.yml"))
        for path in ('"docs/site/**"', '"scripts/site/**"', '"docs/index-calculus-scoreboard.html"'):
            # Once for pull_request, once for push.
            self.assertEqual(workflow.count(path), 2, path)


if __name__ == "__main__":
    unittest.main()
