// Vector-to-PNG renderer for the benchmark browsing copies. Inputs are
// explicit, and the SVG beside each run remains the editable source.
#set page(
  width: float(sys.inputs.at("width")) * 1pt,
  height: float(sys.inputs.at("height")) * 1pt,
  margin: 0pt,
)
#image(sys.inputs.at("image"), width: 100%, height: 100%)
