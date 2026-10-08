#!/usr/bin/env python3
"""Optional browser checks: run with a Python environment containing Playwright."""
from pathlib import Path
import argparse
import tempfile


def main():
    from playwright.sync_api import sync_playwright, expect
    from build import build

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screenshots', type=Path)
    args = parser.parse_args()
    if args.screenshots:
        args.screenshots.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix='crypto-algorithm-lab-') as tmp:
        site = Path(tmp) / 'site'
        build(str(site), str(root))
        url = (site / 'scoreboard/algorithm-lab.html').as_uri()
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page(viewport={'width': 1440, 'height': 1000})
            errors = []
            requests = []
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.on('request', lambda request: requests.append(request.url))
            page.goto(url)
            expect(page.locator('#decompositions')).to_contain_text('F1 + F2')
            page.locator('#matrix-finish').click()
            expect(page.locator('#matrix-solution')).to_contain_text('ℓ1 = 7 · ℓ2 = 4 · ℓ3 = 1 · ℓ4 = 2')
            expect(page.locator('#matrix-solution')).to_contain_text('[11]G = (13, 10): verified')
            for size in ('2', '4', '9'):
                page.locator('#base-size').select_option(size)
                expect(page.locator('.point-chip')).to_have_count(int(size))
                expect(page.locator('#matrix-back')).to_be_disabled()
                page.locator('#matrix-finish').click()
                expect(page.locator('#matrix-solution')).to_contain_text('Every recovered logarithm passes')
            page.locator('#base-size').select_option('2')
            page.locator('#target-k').evaluate("el => { el.value = '1'; el.dispatchEvent(new Event('input', {bubbles:true})); }")
            expect(page.locator('#decompositions')).to_contain_text('No two-point decomposition')
            page.locator('#find-hit').click()
            expect(page.locator('#decompositions')).to_contain_text('verified decomposition')
            page.locator('#matrix-case').select_option('dependent')
            page.locator('#matrix-finish').click()
            expect(page.locator('#matrix-solution')).to_contain_text('Underdetermined')
            page.locator('#matrix-case').select_option('inconsistent')
            page.locator('#matrix-finish').click()
            expect(page.locator('#matrix-solution')).to_contain_text('Rejected')
            page.locator('#matrix-back').click()
            expect(page.locator('#matrix-next')).to_be_enabled()
            page.locator('#matrix-reset').click()
            expect(page.locator('#matrix-back')).to_be_disabled()
            page.locator('#base-size').select_option('4')
            page.locator('#matrix-case').select_option('full')
            page.locator('#target-k').evaluate("el => { el.value = '11'; el.dispatchEvent(new Event('input', {bubbles:true})); }")
            page.locator('#matrix-next').click()
            if args.screenshots:
                page.locator('#factor-base').screenshot(path=str(args.screenshots / 'factor-base-desktop.png'))
                page.locator('#linear-algebra').screenshot(path=str(args.screenshots / 'linear-algebra-desktop.png'))
            # Observe a real backtrack before running the full combination set.
            found_backtrack = False
            for _ in range(30):
                if page.locator('#sat-event').inner_text() == 'Backtrack':
                    found_backtrack = True
                    if args.screenshots:
                        page.locator('#sat-solver').screenshot(path=str(args.screenshots / 'sat-backtrack-desktop.png'))
                    break
                if page.locator('#sat-next').is_disabled():
                    break
                page.locator('#sat-next').click()
            assert found_backtrack
            for kind in ('gates', 'search', 'unsat'):
                page.locator('#sat-case').select_option(kind)
                for encoding in ('native', 'cnf'):
                    page.locator('#sat-encoding').select_option(encoding)
                    for first in ('0', '1'):
                        page.locator('#sat-first').select_option(first)
                        page.locator('#sat-finish').click()
                        expect(page.locator('#sat-event')).to_have_text('UNSAT — no model' if kind == 'unsat' else 'SAT — model found')
                        expect(page.locator('#truth-table tbody tr')).to_have_count(16)
                        expect(page.locator('#sat-next')).to_be_disabled()
            with page.expect_download() as download_info:
                page.locator('#sat-download').click()
            download = download_info.value
            assert download.suggested_filename == 'algorithm-lab-unsat.cnf'
            assert 'p cnf 4 ' in Path(download.path()).read_text()
            page.locator('#sat-reset').click()
            expect(page.locator('#sat-back')).to_be_disabled()
            for width in (390, 320):
                page.set_viewport_size({'width': width, 'height': 844})
                page.locator('#base-size').select_option('9')
                assert page.evaluate('document.documentElement.scrollWidth <= window.innerWidth'), f'Page overflow at {width}px'
                page.locator('#matrix-finish').click()
                expect(page.locator('#matrix-solution')).to_contain_text('Every recovered logarithm passes')
                page.locator('#sat-case').select_option('gates')
                page.locator('#sat-finish').click()
                expect(page.locator('#sat-event')).to_have_text('SAT — model found')
                if args.screenshots and width == 390:
                    page.screenshot(path=str(args.screenshots / 'algorithm-lab-mobile.png'), full_page=True)
            # Built performance page links to the correct lab route.
            page.goto((site / 'scoreboard/performance-gains.html').as_uri())
            page.get_by_role('link', name='Explore factor bases', exact=True).click()
            assert page.url == url + '#factor-base'
            no_js = browser.new_context(java_script_enabled=False)
            fallback = no_js.new_page()
            fallback.goto(url)
            expect(fallback.locator('noscript p.notice')).to_be_visible()
            fallback.get_by_text('How the example chooses and checks points', exact=True).click()
            expect(fallback.get_by_text('Worked example at B = 4:', exact=False)).to_be_visible()
            assert not errors, errors
            assert not any(request.startswith(('http:', 'https:')) for request in requests), 'Demo made a network request'
            browser.close()
    print('PASS: desktop/mobile UI, linked matrix state, backtracking, 12 SAT configurations, download, Pages navigation, no-JS fallback, no network requests.')


if __name__ == '__main__':
    main()
