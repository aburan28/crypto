#!/bin/sh
# Re-render the journey graphs and PROGRAM.pdf from recorded entries (layout only).
# Needs python3 with matplotlib and typst (pip install typst matplotlib in a scratch venv).
set -e
here=$(cd "$(dirname "$0")" && pwd)
PY=${PY:-python3}
"$PY" -I "$here/journey/render_journey.py" "$here"
"$PY" -I -c "import typst,sys; typst.compile(sys.argv[1]+'/PROGRAM.typ', output=sys.argv[1]+'/PROGRAM.pdf')" "$here"
