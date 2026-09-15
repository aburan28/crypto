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
            "assets/rho-gpu.js",
            "assets/rho-gpu-worker.js",
            "assets/rho-gpu-host.js",
            "assets/rho-gpu.wgsl",
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

    def test_rho_engine_imports_and_fetches_resolve_beside_the_worker(self):
        # The worker is a module that imports the host arithmetic and fetches
        # the shader by relative URL, both resolved against its own location.
        # Publishing it without one of them gives a page whose GPU control
        # fails only when someone switches it on.
        worker = read(os.path.join(self.out, "assets", "rho-gpu-worker.js"))
        targets = re.findall(r'from "\./([^"]+)"', worker) + re.findall(r'fetch\("\./([^"]+)"', worker)
        self.assertIn("rho-gpu-host.js", targets)
        self.assertIn("rho-gpu.wgsl", targets)
        for target in targets:
            self.assertTrue(
                os.path.exists(os.path.join(self.out, "assets", target)),
                "the rho worker loads %s, which the build does not publish beside it" % target,
            )

    def test_rho_engine_is_off_until_the_visitor_switches_it_on(self):
        # Spending a visitor's GPU without being asked is the one thing this
        # feature must never do, and it is one attribute away at all times:
        # a `checked` on the control, or a start that does not read the stored
        # choice. Pin both, and the default-off of the background-tab option.
        page = read(os.path.join(self.out, "index.html"))
        controls = re.findall(r'<input type="checkbox" id="rho-(?:run|bg)"[^>]*>', page)
        self.assertEqual(len(controls), 2, controls)
        for control in controls:
            self.assertNotIn("checked", control, control)
        script = read(os.path.join(self.out, "assets", "rho-gpu.js"))
        self.assertIn("toggle.checked = !!prefs.on;", script)
        self.assertIn("bg.checked = !!prefs.bg;", script)
        # Hidden tabs pause unless the visitor asked otherwise.
        self.assertIn("return document.hidden && !bg.checked;", script)

    def test_rho_engine_reports_no_rate_before_its_device_is_verified(self):
        # A measurement from a device that does not agree with the host
        # reference is not a measurement. The worker replays its first device
        # steps in BigInt and throws before it posts a "ready".
        worker = read(os.path.join(self.out, "assets", "rho-gpu-worker.js"))
        self.assertIn("self-test failed: device and host disagree", worker)
        self.assertLess(
            worker.index("self-test failed: device and host disagree"),
            worker.index('type: "ready"'),
            "the self-test must run before the engine reports itself ready",
        )

    def test_pages_do_not_render_fetched_data_through_innerhtml(self):
        # The dashboard builds rows with textContent because worker_id is
        # self-reported by the walkers; keep it that way.
        for rel in ("index.html", "status/index.html"):
            page = read(os.path.join(self.out, rel))
            self.assertNotIn("innerHTML", page, rel)

    def test_landing_page_operations_total_cites_the_campaign_constants(self):
        # The landing page turns the distinguished-point count into a
        # group-operation total, so it depends on two figures it does not
        # measure: the campaign's report rate at HW(x) <= 34, and the
        # expected cost of a collision. Both are documented in ecc2k130/;
        # pin them together so the page cannot drift from what it cites.
        page = read(os.path.join(self.out, "index.html"))
        campaign = read(os.path.join(ROOT, "ecc2k130", "aws", "README.md"))
        self.assertIn('id="live-ops"', page)
        for exponent in ("2^25.27", "2^60.9"):
            self.assertIn(exponent, page, exponent)
            self.assertIn(exponent, campaign, exponent)

    def test_dashboard_iteration_total_matches_the_landing_page(self):
        # Both published pages turn the same point count into a walk total,
        # so a constant edited on one of them only would publish two different
        # iteration counts for one campaign. Pin the two pages to each other
        # and to the campaign document both of them cite.
        dashboard = read(os.path.join(self.out, "status", "index.html"))
        landing = read(os.path.join(self.out, "index.html"))
        campaign = read(os.path.join(ROOT, "ecc2k130", "aws", "README.md"))
        self.assertIn('id="ops-value"', dashboard)
        for exponent in ("2^25.27", "2^60.9"):
            for name, page in (("dashboard", dashboard), ("landing", landing), ("campaign", campaign)):
                self.assertIn(exponent, page, "%s is missing %s" % (name, exponent))
        for name, page in (("dashboard", dashboard), ("landing", landing)):
            self.assertIn("ITERATIONS_PER_DP_LOG2 = 25.27;", page, name)
            self.assertIn('CAMPAIGN = "ecc2k-130";', page, name)
        self.assertIn("EXPECTED_ITERATIONS_LOG2 = 60.9;", dashboard)

    def test_both_pages_read_the_measured_rate_through_one_shared_block(self):
        # Two pages render the same snapshot's rate, so a fix applied to one
        # copy and not the other publishes two different walk rates for one
        # campaign. The block is small enough to mirror and too important to
        # let drift, so pin the copies byte for byte.
        start = "  // --- mirrored block:"
        end = "  // --- end mirrored block ---"

        def mirrored(rel):
            page = read(os.path.join(self.out, rel))
            self.assertIn(start, page, rel)
            self.assertIn(end, page, rel)
            return page[page.index(start):page.index(end) + len(end)]

        dashboard = mirrored("status/index.html")
        landing = mirrored("index.html")
        self.assertEqual(dashboard, landing)
        for name in ("function reportedIterations", "function measuredRate", "function formatRate"):
            self.assertIn(name, dashboard, name)
        # B it/s is the unit the campaign quotes a GPU in (ecc2k130/aws/README.md).
        self.assertIn('" B it/s"', dashboard)

    def test_pages_prefer_the_counted_iteration_total_over_the_derived_one(self):
        # The derived total is the point count times the interval for
        # HW(x) <= 34, and this campaign distinguishes at HW(x) <= 32, so the
        # derivation reads about 2^2.7 low. Both pages must take the walkers'
        # own count when the snapshot carries one, and both must still be able
        # to fall back for a snapshot that does not.
        # The dashboard names the interval the running cutoff measures at, so
        # that figure has to exist in the campaign document it comes from.
        campaign = read(os.path.join(ROOT, "ecc2k130", "aws", "README.md"))
        dashboard = read(os.path.join(self.out, "status", "index.html"))
        for cited in ("HW(x) &le; 32", "2^27.9"):
            self.assertIn(cited, dashboard, cited)
        self.assertIn("2^27.9", campaign)
        for rel in ("status/index.html", "index.html"):
            page = read(os.path.join(self.out, rel))
            self.assertIn("var reported = reportedIterations(status);", page, rel)
            body = page[page.index("var reported = reportedIterations(status);"):]
            body = body[:body.index("\n  }")]
            self.assertLess(
                body.index("return"),
                body.index("ITERATIONS_PER_DP_LOG2"),
                "%s applies the fallback interval before checking for a counted total" % rel,
            )

    def test_dashboard_publishes_the_rate_and_the_span_it_measured(self):
        # A rate without its window is not checkable: 15 minutes of
        # checkpoints and an hour of them are different measurements, and the
        # short one swings by a third on this fleet. The value, the span and
        # the iteration series behind them all have to be on the page.
        page = read(os.path.join(self.out, "status", "index.html"))
        self.assertIn('id="walk-rate"', page)
        self.assertIn('id="walk-rate-sub"', page)
        self.assertIn("function drawRate", page)
        self.assertIn("drawRate(status)", page)
        self.assertIn("status.walk_rate.window_seconds", page)
        self.assertIn("mean over the last ", page)
        # The history table carries the totals the rate is the difference of.
        self.assertIn("num(point.iterations)", page)
        self.assertIn('<th scope="col" class="num">iterations</th>', page)

    def test_dashboard_eta_runs_on_the_measured_rate_and_falls_back_to_dps(self):
        # The ETA is the remaining expected work divided by a rate, and which
        # rate it is matters: the measured one when two snapshots carry an
        # iteration total, and otherwise the last-hour DP amount at 2^25.27
        # per point. Pin the formula and the order, not a rendered duration.
        page = read(os.path.join(self.out, "status", "index.html"))
        self.assertIn('id="eta-value"', page)
        self.assertIn("function opsPerSecond", page)
        self.assertIn("function etaSeconds", page)
        self.assertIn("function drawEta", page)
        ops = page[page.index("function opsPerSecond"):page.index("function etaSeconds")]
        self.assertIn("var measured = measuredRate(status);", ops)
        self.assertIn("dpsLastHour * Math.pow(2, ITERATIONS_PER_DP_LOG2)) / 3600", ops)
        self.assertLess(
            ops.index("measuredRate(status)"),
            ops.index("dps_last_hour"),
            "the ETA prefers the DP derivation over the measured rate",
        )
        # And never mixes them: a counted total over a derived rate is an ETA
        # six times too long, which is worse than no ETA.
        self.assertIn("if (reportedIterations(status) !== null) return null;", ops)
        self.assertIn(
            "Math.pow(2, EXPECTED_ITERATIONS_LOG2) - Math.pow(2, log2ops)",
            page,
        )
        self.assertIn("drawEta(status)", page)
        # Without a measured rate the last-hour card foot still surfaces the
        # ops/h that amount implies, so a change in the interval's count is
        # visible as ops/h rather than only as a point count.
        self.assertIn("ops/h at interval 2^", page)

    def test_dashboard_progress_bar_is_linear_in_work(self):
        # The share of 2^60.9 walked so far is around 2^-21, so a bar drawn
        # from the ratio of the EXPONENTS would read about two thirds full
        # while the campaign has done a millionth of a millionth of the work.
        # The visible progress bar must therefore be filled from the ratio of
        # the work itself, and the log-scale bar beside it must say in the
        # page that it is not progress. Both are easy to "fix" into a lie by
        # someone making the bar look better, so pin them here.
        page = read(os.path.join(self.out, "status", "index.html"))
        self.assertIn("var share = Math.pow(2, log2ops - EXPECTED_ITERATIONS_LOG2);", page)
        self.assertIn("var percent = Math.min(100, share * 100);", page)
        self.assertIn('fill.style.width = percent + "%";', page)
        self.assertIn("this bar is not progress", page)
        # A minimum width on the fill would draw a share that is not there.
        match = re.search(r"\.bar-fill \{(.*?)\}", read(os.path.join(self.out, "status", "style.css")), re.S)
        self.assertIsNotNone(match, "no .bar-fill rule in the dashboard stylesheet")
        self.assertNotIn("min-width", match.group(1))

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

    def test_research_evidence_is_linked_absolutely(self):
        # research/ is not part of the site artifact, so a relative
        # ../research/<path> link resolves in the working tree and 404s once
        # published under /scoreboard/. It has broken publication twice; the
        # generic link check above catches it, but only by filename, so this
        # names the cause and the fix. See scripts/site/README.md.
        relative = []
        for rel in ("index.html", "404.html", "status/index.html", "scoreboard/index.html"):
            page = read(os.path.join(self.out, rel))
            for href in re.findall(r'(?:href|src)="([^"]+)"', page):
                if re.match(r"(?:\.\./)*research/", href):
                    relative.append("%s -> %s" % (rel, href))
        self.assertEqual(
            relative,
            [],
            "link research evidence as https://github.com/aburan28/crypto/blob/main/research/<path>, "
            "not relatively: %s" % relative,
        )

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
