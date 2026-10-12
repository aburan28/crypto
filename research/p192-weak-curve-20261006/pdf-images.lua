-- Keep SVG as the canonical Markdown rendering while selecting vector PDF
-- siblings for the XeLaTeX build. The siblings are generated from the same
-- editable sources immediately before REPORT.pdf.
function Image(image)
  if FORMAT:match("latex") then
    image.src = image.src:gsub("%.svg$", ".pdf")
  end
  return image
end
