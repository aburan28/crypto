"""Document rendering only; all research values come from native analysis."""
from pathlib import Path
import io
import shutil

import cairosvg
from PIL import Image
import typst

STUDY = Path(__file__).resolve().parent
REPO = STUDY.parents[1]
GALLERY = REPO / "benchmarks/20261014-conductor-gap-instructions"
DOCUMENTS = [
    (STUDY, "report", "instruction-ratios"),
    (STUDY / "reports/ell8209-q29", "ell8209-q29", "instructions"),
    (STUDY / "reports/ell131101-q37", "ell131101-q37", "instructions"),
    (STUDY / "attempt-unpreregistered", "attempt-unpreregistered", "attempt"),
    (STUDY / "attempt-phase-dumps", "attempt-phase-dumps", "attempt"),
]


def main():
    GALLERY.mkdir(parents=True, exist_ok=True)
    for directory, name, graph in DOCUMENTS:
        typst.compile(str(directory / "report.typ"), output=str(directory / "report.pdf"))
        shutil.copyfile(directory / "report.pdf", GALLERY / f"{name}.pdf")
        shutil.copyfile(directory / f"{graph}.svg", GALLERY / f"{name}.svg")
        png = cairosvg.svg2png(url=str(directory / f"{graph}.svg"), scale=150 / 96)
        with Image.open(io.BytesIO(png)) as image:
            image.save(GALLERY / f"{name}.png", dpi=(150, 150))
        print(f"rendered {directory.relative_to(REPO)}")


if __name__ == "__main__":
    main()
