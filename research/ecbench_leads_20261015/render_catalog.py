"""Render existing catalog prose/HTML; new experimental figures come from Rust.

Do not call build(): it derives the historical research data. Only page and
markdown lay out the already frozen native-generated leaderboard document.
"""
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
import build_ic_leaderboard as display

doc=json.loads(display.OUT_JSON.read_text())
display.OUT_HTML.write_text(display.page(doc,standalone=True))
display.OUT_MD.write_text(display.markdown(doc))
