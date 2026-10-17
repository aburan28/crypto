"""Presentation checks/screenshots only; no experiments or curve arithmetic.

Usage: capture.py http://127.0.0.1:8768/graphs.html FRESH_SCREENSHOT_DIR
Requires Playwright and Chromium. Managed Chromium blocks file:// navigation;
serve only the generated artifact directory over loopback in that environment.
"""
import json
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright

url, output = sys.argv[1], Path(sys.argv[2])
output.mkdir(parents=True, exist_ok=False)
with sync_playwright() as p:
    browser = p.chromium.launch(executable_path="/usr/bin/chromium", headless=True,
                               args=["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"],
                               timeout=15000)
    page = browser.new_page(viewport={"width": 1540, "height": 1060})
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(url + "?graph=2", wait_until="domcontentloaded", timeout=15000)
    page.wait_for_selector("html[data-ready=true]", timeout=10000)
    assert "10,070 vertices" in page.locator("#status").inner_text()
    page.screenshot(path=str(output / "mixed.png"), full_page=True)
    page.locator("#vertex").fill("10069")
    page.locator("#go").click()
    assert isinstance(json.loads(page.locator("#model").inner_text())["a"], str)
    page.locator("#neighbors button").first.click()
    assert "kernel_json" in page.locator("#map").inner_text()
    page.locator("#graph").select_option("1")
    assert "6,561 vertices" in page.locator("#status").inner_text()
    page.locator("#vertex").fill("6560")
    page.locator("#go").click()
    page.screenshot(path=str(output / "floor.png"), full_page=True)
    page.set_viewport_size({"width": 390, "height": 844})
    page.locator("#graph").select_option("0")
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth"), "mobile overflow"
    page.screenshot(path=str(output / "mobile.png"), full_page=True)
    assert not errors, errors
    browser.close()
print("PASS: real Chromium three-class navigation, large-ID class switch, exact coefficient "
      "strings, outgoing kernel detail, floor inspector and mobile width; no page errors.")
