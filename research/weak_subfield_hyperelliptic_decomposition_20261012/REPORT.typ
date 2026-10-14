// Results report: renders README.md (the source of truth) with its figures.
#import "@preview/cmarker:0.1.6"
#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 7.5pt)
#show raw.where(block: false): set text(size: 0.92em)
#cmarker.render(read("README.md"), scope: (image: (path, alt: none) => image(path, width: 100%)))
