-- Keep SVG as the canonical Markdown rendering while selecting the vector
-- PDF sibling, generated from the same DOT source, for the XeLaTeX build.
function Image(image)
  if FORMAT:match("latex") then
    image.src = image.src:gsub("%.svg$", ".pdf")
  end
  return image
end
