// Build: pandoc -f gfm -t typst REPORT.md -o report.typ; typst compile report-main.typ REPORT.pdf
#set page(paper: "a4", margin: (x: 1.4cm, y: 1.6cm), numbering: "1")
#set text(size: 8.5pt, font: ("Libertinus Serif", "New Computer Modern"))
#show raw: set text(size: 7.5pt)
#show table: set text(size: 7pt)
#show image: it => block(width: 100%, it)
#set table(inset: 3pt)
