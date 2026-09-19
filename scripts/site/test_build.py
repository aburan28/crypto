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

FOREST_DIR = os.path.join(ROOT, "docs", "ecc2k130-status", "walk-forest")


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
            "favicon.svg",
            "assets/site.css",
            "assets/rho-gpu.js",
            "assets/rho-gpu-worker.js",
            "assets/rho-gpu-host.js",
            "assets/rho-gpu.wgsl",
            "scoreboard/index.html",
            "scoreboard/algorithm-lab.html",
            "scoreboard/algorithm-lab/core.js",
            "scoreboard/algorithm-lab/ui.js",
            "scoreboard/algorithm-lab/style.css",
            "scoreboard/algorithm-lab/README.md",
            "scoreboard/performance-gains.html",
            "scoreboard/performance-gains/summary.png",
            "scoreboard/performance-gains/summary.pdf",
            "scoreboard/performance-gains/summary.svg",
            "scoreboard/performance-gains/data.json",
            "scoreboard/performance-gains/comparisons.csv",
            "status/index.html",
            "status/style.css",
            "status/walk-forest.svg",
            "status/walk-forest.json",
            "status/walk-forest-gf2-23.json",
            "status/walk-forest.js",
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
        for name in (
            "function reportedIterations",
            "function measuredRate",
            "function walkingSlots",
            "function formatRate",
            "function campaignDisplayState",
        ):
            self.assertIn(name, dashboard, name)
        # B it/s is the unit the campaign quotes a GPU in (ecc2k130/aws/README.md).
        self.assertIn('" B it/s"', dashboard)
        self.assertIn("status.work.walking_slots", dashboard)

    def test_pages_show_walking_slots_as_gpus_running(self):
        # status.workers is lifetime DISTINCT worker_id from the DP table;
        # GPUs that are actually walking are work.walking_slots from the
        # checkpoint feed. Publishing the lifetime count in the headline made
        # a three-GPU fleet read as three thousand.
        dashboard = read(os.path.join(self.out, "status", "index.html"))
        landing = read(os.path.join(self.out, "index.html"))
        self.assertIn('id="gpus"', dashboard)
        self.assertIn('id="gpus-foot"', dashboard)
        self.assertIn(">GPUs running<", dashboard)
        self.assertIn("function drawGpus", dashboard)
        self.assertIn("drawGpus(status)", dashboard)
        self.assertNotIn('el("gpus").textContent = num(status.workers)', dashboard)
        self.assertIn('id="live-gpus"', landing)
        self.assertIn(">GPUs running<", landing)
        self.assertIn("walkingSlots(status)", landing)
        self.assertIn("campaignDisplayState(status)", landing)
        self.assertNotIn("live-workers", landing)
        # Lifetime contributors stay in the workers table, not the GPU card.
        self.assertIn('id="workers-note"', dashboard)
        self.assertIn("Lifetime distinguished-point contributors", dashboard)

    def test_pages_prefer_the_counted_iteration_total_over_the_derived_one(self):
        # The derived total is the point count times the interval for
        # HW(x) <= 34, and this campaign distinguishes at HW(x) <= 32, so the
        # derivation reads about 2^3.1 low. Both pages must take the walkers'
        # own count when the snapshot carries one, and both must still be able
        # to fall back for a snapshot that does not.
        # The dashboard names the interval the running cutoff measures at, so
        # that figure has to exist in the campaign document it comes from.
        campaign = read(os.path.join(ROOT, "ecc2k130", "aws", "README.md"))
        dashboard = read(os.path.join(self.out, "status", "index.html"))
        for cited in ("HW(x) &le; 32", "2^28.4"):
            self.assertIn(cited, dashboard, cited)
        self.assertIn("2^28.4", campaign)
        # The retired 2026-09-15 ratio must not come back as the live figure.
        self.assertNotIn("2^27.9</span> iterations per point", dashboard)
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

    def test_dashboard_coordinate_table_cites_the_priced_note(self):
        # The coordinate-system table prices one update per point
        # representation. Every figure on it comes from
        # ecc2k130/LAMBDA-PROJECTIVE.md and nowhere else (AGENTS.md section 7:
        # the page cites, never computes), so pin each carry-less figure on
        # the page to that note, the page's link to the note, and the landing
        # page's entry pointing at both.
        dashboard = read(os.path.join(self.out, "status", "index.html"))
        landing = read(os.path.join(self.out, "index.html"))
        # the note bolds its reference row, so drop the markers before matching
        note = read(os.path.join(ROOT, "ecc2k130", "LAMBDA-PROJECTIVE.md")).replace("*", "")
        self.assertIn('id="h-coords"', dashboard)
        section = dashboard[dashboard.index('id="h-coords"'):dashboard.index('id="h-workers"')]
        figures = re.findall(r'<td class="num clmad">([^<]+)</td>', section)
        self.assertGreaterEqual(len(figures), 8, figures)
        for figure in figures:
            self.assertIn("| %s |" % figure, note, "dashboard prices %s clmad, which the note does not carry" % figure)
        self.assertIn("ecc2k130/LAMBDA-PROJECTIVE.md", section)
        self.assertIn("ecc2k130/LAMBDA-PROJECTIVE.md", landing)
        self.assertIn('href="./status/#h-coords"', landing)

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
        # nine times too long, which is worse than no ETA.
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

    def test_dashboard_shows_the_work_share_as_a_number_not_a_bar(self):
        # The share of 2^60.9 walked was 2^-21 when the page was first
        # published: a bar filled from it is empty, and a bar filled from the
        # ratio of the exponents reads two thirds full. Both are lies waiting
        # to be drawn, so the share is a number and the page has no progress
        # bar to make look better.
        page = read(os.path.join(self.out, "status", "index.html"))
        self.assertIn("var share = Math.pow(2, log2ops - EXPECTED_ITERATIONS_LOG2);", page)
        self.assertIn('el("ops-percent").textContent = percentText(share);', page)
        self.assertNotIn('role="progressbar"', page)
        self.assertNotIn("ops-fill", page)

    def test_dashboard_odds_panel_is_the_birthday_bound_with_the_quoted_mean(self):
        # What rho accumulates is the chance a collision has already happened.
        # The page draws P(W) = 1 - exp(-pi W^2 / 4 E^2): the birthday bound
        # parameterised so its MEAN is E = 2^60.9, the expectation the page
        # already quotes, so the curve and the ETA cannot disagree. Pin the
        # formula, its inverse, the conditional month tile, the doubled-fleet
        # curve and the "not a deadline" wording; a rounder-looking curve is
        # one edit away.
        page = read(os.path.join(self.out, "status", "index.html"))
        self.assertIn("function collisionOdds", page)
        self.assertIn("return -Math.expm1(-Math.PI / 4 * ratio * ratio);", page)
        self.assertIn("function workForOdds", page)
        self.assertIn("Math.sqrt(-4 * Math.log(1 - p) / Math.PI)", page)
        self.assertIn("function drawOdds", page)
        self.assertIn("drawOdds(status)", page)
        self.assertIn("(oddsMonth - odds0) / (1 - odds0)", page)
        self.assertIn('shape("path", { d: curve(2), "class": "line alt" })', page)
        self.assertIn("The median is not a deadline", page)
        # Axis dates and guide marks no longer share one strip under the plot —
        # that stacking is what clipped the labels on a phone-width viewport.
        self.assertIn("PAD_B = 52", page)
        self.assertIn('label(gx, H - 22, xLabel(work, mark), "middle", "tick")', page)
        self.assertIn('label(gx + 4, y(p) - 6, mark, "start", "tick")', page)
        css = read(os.path.join(self.out, "status", "style.css"))
        self.assertIn(".odds-chart svg", css)
        self.assertIn("min-width: 0", css)
        # Walk forest fits the column on a phone. A 640px min-width made the
        # caption lay out at that width, so every line clipped; the SVG
        # thickens its own strokes below 720px instead.
        self.assertIn(".figure img", css)
        self.assertIn("min-width: 0", css)
        self.assertNotIn("min-width: 640px", css)
        self.assertNotIn("Swipe to pan the forest", css)
        self.assertIn("figure-frame", page)
        self.assertIn("@media(max-width:720px)", read(os.path.join(self.out, "status", "walk-forest.svg")))
        for ident in ("odds-now", "odds-month", "odds-median", "odds-ninety", "odds-chart"):
            self.assertIn('id="%s"' % ident, page, ident)
        # The same law in Python: mean E, median 0.94 E, ninety at 1.71 E.
        import math
        odds = lambda ratio: -math.expm1(-math.pi / 4 * ratio * ratio)
        work_for = lambda p: math.sqrt(-4 * math.log(1 - p) / math.pi)
        self.assertAlmostEqual(odds(work_for(0.5)), 0.5)
        self.assertAlmostEqual(work_for(0.5), 0.939, places=3)
        self.assertAlmostEqual(work_for(0.9), 1.712, places=3)
        mean = sum((1 - odds(k / 1000.0)) * 0.001 for k in range(6000))
        self.assertAlmostEqual(mean, 1.0, places=2)
        # The page's copy of the two quantiles agrees.
        self.assertIn("0.94 E", page)
        self.assertIn("1.71 E", page)

    # ---- the walk-forest figure ------------------------------------------
    # docs/ecc2k130-status/walk-forest.svg is generated from the sampled,
    # hashed trails committed beside it: real ECC2K-130 walks the client
    # made, replayed with its own kernel.  These hold the published figure,
    # the page's caption, the endpoints the trails must land on and the
    # scalar reference's replay to one another.

    def _forest(self):
        import walk_forest
        forest = walk_forest.read_trails(os.path.join(FOREST_DIR, "trails.txt"))
        corpus = walk_forest.read_corpus(os.path.join(FOREST_DIR, "forest.hashes"))
        return walk_forest, forest, corpus

    def test_walk_forest_is_sampled_from_the_challenge_curve(self):
        walk_forest, forest, corpus = self._forest()
        p = forest.params
        self.assertEqual(p["curve"], "131")
        self.assertEqual(p["mode"], "sample")
        self.assertEqual(p["source"], "corpus")
        self.assertEqual(p["hash"], "sha256-16")
        # Every seed the sampler was given came from a client record, and
        # every trail it kept ended on the orbit its record names.
        self.assertEqual(int(p["checked"]), int(p["drawn"]))
        self.assertEqual(int(p["drawn"]), len(forest.walks))
        self.assertLessEqual(int(p["drawn"]), int(p["tried"]))
        for steps in forest.steps:
            self.assertLessEqual(steps, int(p["cap"]))
        # Orbits are named, never shown: a 16-hex-digit name per node.
        for seed, orbits in forest.walks:
            for orbit in orbits:
                self.assertRegex(orbit, r"^[0-9a-f]{16}$")

    def test_walk_forest_trails_end_on_their_recorded_endpoints(self):
        walk_forest, forest, corpus = self._forest()
        self.assertEqual(walk_forest.bind_to_corpus(forest, corpus), len(forest.walks))
        self.assertEqual(len(corpus), len(forest.walks))

    def test_walk_forest_kernel_replay_agrees_with_the_scalar_reference(self):
        # reference-check.txt is the scalar reference's replay of the first
        # records, hashed and sampled the same way (walk_forest.py
        # --hash-trails).  The kernel's trails for those walks must match it
        # orbit for orbit.
        walk_forest, forest, corpus = self._forest()
        reference = walk_forest.read_trails(os.path.join(FOREST_DIR, "reference-check.txt"))
        self.assertGreater(len(reference.walks), 0)
        self.assertEqual(reference.every, forest.every)
        for (rseed, rorbits), (seed, orbits), rsteps, steps in zip(
                reference.walks, forest.walks, reference.steps, forest.steps):
            self.assertEqual(rseed, seed)
            self.assertEqual(rsteps, steps, seed)
            self.assertEqual(rorbits, orbits, seed)

    def test_walk_forest_figure_is_what_the_trails_render_to(self):
        walk_forest, forest, corpus = self._forest()
        svg, _ = walk_forest.build(os.path.join(FOREST_DIR, "trails.txt"), os.path.join(FOREST_DIR, "forest.hashes"))
        published = read(os.path.join(self.out, "status", "walk-forest.svg"))
        self.assertEqual(published, svg, "walk-forest.svg is stale: regenerate it (see walk-forest/README.md)")
        # And the drawing carries exactly the forest: one circle per drawn
        # orbit, one filled circle per distinguished point, one edge per
        # stride.
        self.assertEqual(published.count("<circle "), len(forest.order))
        self.assertEqual(published.count('class="dp"'), len(forest.roots))
        plain = published.split('<g class="e">', 1)[1].split("</g>", 1)[0]
        self.assertEqual(plain.count("M"), len(forest.succ))
        self.assertIn(forest.header, published)

    def test_walk_forest_caption_states_the_counts_it_draws(self):
        walk_forest, forest, corpus = self._forest()
        page = read(os.path.join(self.out, "status", "index.html"))
        meetings = sum(1 for preds in forest.pred.values() if len(preds) > 1)
        expected = {
            "wf-walks": len(forest.walks),
            "wf-tried": int(forest.params["tried"]),
            "wf-cap": int(forest.params["cap"]),
            "wf-every": forest.every,
            "wf-orbits": len(forest.order),
            "wf-steps": walk_forest.total_steps(forest),
            "wf-longest": max(forest.steps),
            "wf-dps": len(forest.roots),
            "wf-meetings": meetings,
            "wf-weight": int(forest.params["dp-weight"]),
        }
        for ident, value in expected.items():
            match = re.search(r'id="%s">([^<]+)<' % ident, page)
            self.assertIsNotNone(match, ident)
            self.assertEqual(int(match.group(1).replace(",", "")), value, ident)
        self.assertIn('src="./walk-forest.svg"', page)

    def test_walk_forest_graphs_are_what_the_trails_export_to(self):
        # The explorer draws walk-forest.json (the real curve) and
        # walk-forest-gf2-23.json (the test curve) with the static figure's
        # own layout, so each must be byte-identical to what its trails
        # export to, carry every node, edge and walk, and name nodes the way
        # the trails do: hash prefixes on the challenge curve.
        import json
        import walk_forest
        sets = (
            ("walk-forest.json", "trails.txt", "forest.hashes", True),
            ("walk-forest-gf2-23.json", os.path.join("gf2-23", "trails.txt"), os.path.join("gf2-23", "forest.bin"), False),
        )
        for published, trails, corpus, hashed in sets:
            text, forest = walk_forest.build_graph(
                os.path.join(FOREST_DIR, trails), os.path.join(FOREST_DIR, corpus),
                title=json.loads(read(os.path.join(self.out, "status", published)))["title"])
            self.assertEqual(read(os.path.join(self.out, "status", published)), text,
                             "%s is stale: regenerate it (see walk-forest/README.md)" % published)
            graph = json.loads(text)
            self.assertEqual(len(graph["nodes"]), len(forest.order), published)
            self.assertEqual(len(graph["edges"]), len(forest.succ), published)
            self.assertEqual(len(graph["walks"]), len(forest.walks), published)
            self.assertEqual(sum(1 for n in graph["nodes"] if n[3]), len(forest.roots), published)
            self.assertEqual(graph["counts"]["meetings"], sum(1 for p in forest.pred.values() if len(p) > 1), published)
            if hashed:
                for node in graph["nodes"]:
                    self.assertRegex(node[0], r"^[0-9a-f]{16}$")
            for walk in graph["walks"]:
                for a, b in zip(walk["nodes"], walk["nodes"][1:]):
                    self.assertEqual(forest.succ[forest.order[a]], forest.order[b], published)

    def test_walk_forest_explorer_is_wired_and_degrades_to_the_figure(self):
        # The script fetches the two graphs by relative URL beside the page,
        # replaces the image only once a graph has drawn, and never renders
        # fetched data through innerHTML (node names are data, hashed or not).
        page = read(os.path.join(self.out, "status", "index.html"))
        script = read(os.path.join(self.out, "status", "walk-forest.js"))
        self.assertIn('<script src="./walk-forest.js" defer></script>', page)
        for ident in ("forest", "forest-dataset", "forest-play", "forest-reset", "forest-readout"):
            self.assertIn('id="%s"' % ident, page, ident)
        self.assertIn('<img src="./walk-forest.svg"', page)
        for target in re.findall(r'file: "\./([^"]+)"', script):
            self.assertTrue(os.path.exists(os.path.join(self.out, "status", target)), target)
        self.assertIn("host.removeChild(image)", script)
        self.assertLess(script.index("index(g);"), script.index("host.removeChild(image)"))
        self.assertNotIn("innerHTML", script)

    # ---- the tool itself, on a curve small enough to check in the clear ---
    # walk-forest/gf2-23/ is the GF(2^23) set: the client's own records for a
    # run, the reference walk of that run's seed schedule, and its endpoints.
    # It exercises the replay and generate modes where orbits can be
    # committed in the clear.

    def test_gf2_23_trails_end_on_their_corpus_records(self):
        import walk_forest
        small = os.path.join(FOREST_DIR, "gf2-23")
        forest = walk_forest.read_trails(os.path.join(small, "trails.txt"))
        corpus = walk_forest.read_corpus(os.path.join(small, "forest.bin"))
        self.assertEqual(forest.params["mode"], "generate")
        self.assertEqual(walk_forest.bind_to_corpus(forest, corpus), len(forest.walks))
        self.assertEqual(len(corpus), len(forest.walks))

    def test_gf2_23_trails_agree_with_the_clients_own_records(self):
        # client-run1.bin is what ecc2k130-cpu itself wrote for the same
        # run-id on this instance before its search solved the instance.
        # Every record of it for a lane the trails cover must name the orbit
        # the trail ends on; the header records how many the generator
        # checked.
        import walk_forest
        small = os.path.join(FOREST_DIR, "gf2-23")
        forest = walk_forest.read_trails(os.path.join(small, "trails.txt"))
        corpus = walk_forest.read_corpus(os.path.join(small, "forest.bin"))
        client = walk_forest.read_corpus(os.path.join(small, "client-run1.bin"))
        overlap = {seed: key for seed, key in client.items() if seed in corpus}
        self.assertGreater(len(overlap), 0, "no client record falls among the covered lanes")
        for seed, key in overlap.items():
            self.assertEqual(walk_forest.canon_hex(key), walk_forest.canon_hex(corpus[seed]), seed)
        self.assertEqual(int(forest.params["checked"]), len(overlap))

    def test_internal_links_resolve(self):
        missing = []
        for rel in ("index.html", "404.html", "status/index.html", "scoreboard/index.html", "scoreboard/performance-gains.html", "scoreboard/algorithm-lab.html"):
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
        for rel in ("index.html", "404.html", "status/index.html", "scoreboard/index.html", "scoreboard/performance-gains.html", "scoreboard/algorithm-lab.html"):
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
