# Chart source: ours / fastest external tool, per prime, from ratios.dat
# (`cargo run --release --example ext_summary -- results/p15-external --dat`).
#   awk [-v YMIN=0.0625 -v DATA=results/p16-external/] -v OP=mul -v XMIN=16 -v XMAX=16384 -v XLABEL="operand length n" -v TITLE="..." -f chart.awk ratios.dat > out.svg
BEGIN { W = 840; H = 480; L = 70; R = 230; T = 46; B = 62; pw = W - L - R; ph = H - T - B
  x0 = log(XMIN); x1 = log(XMAX); y0 = log(YMIN == "" ? 0.25 : YMIN); y1 = log(4)
  col["p61"] = "#6e7781"; col["p127"] = "#1f6feb"; col["p255"] = "#d4580d"; col["p511"] = "#8250df"
  lab["p61"] = "p = 2^61 - 1"; lab["p127"] = "p = 2^127 - 1"; lab["p255"] = "p = 2^255 - 19"; lab["p511"] = "p = CSIDH-512 prime"
  printf "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\" viewBox=\"0 0 %d %d\" font-family=\"sans-serif\" font-size=\"12\">\n<rect width=\"100%%\" height=\"100%%\" fill=\"white\"/>\n", W, H, W, H
  printf "<text x=\"%d\" y=\"24\" font-size=\"15\" font-weight=\"bold\">%s</text>\n", L, TITLE
  n = split(YMIN == 0.0625 ? "0.0625 0.125 0.25 0.5 1 2 4" : "0.25 0.5 1 2 4", ys, " ")
  for (i = 1; i <= n; i++) { y = T + ph * (1 - (log(ys[i]) - y0) / (y1 - y0))
    printf "<line x1=\"%d\" x2=\"%d\" y1=\"%.1f\" y2=\"%.1f\" stroke=\"%s\"/><text x=\"%d\" y=\"%.1f\" text-anchor=\"end\">%s×</text>\n", L, L + pw, y, y, (ys[i] == 1 ? "#555" : "#e5e5e5"), L - 6, y + 4, ys[i] }
  y = T + ph * (1 - (0 - y0) / (y1 - y0))
  printf "<text x=\"%d\" y=\"%.1f\" fill=\"#555\" text-anchor=\"end\">below: ours faster</text><text x=\"%d\" y=\"%.1f\" fill=\"#555\" text-anchor=\"end\">above: external faster</text>\n", L + pw - 4, y + 16, L + pw - 4, y - 6
  for (v = XMIN; v <= XMAX; v *= 4) { x = L + pw * (log(v) - x0) / (x1 - x0)
    printf "<line x1=\"%.1f\" x2=\"%.1f\" y1=\"%d\" y2=\"%d\" stroke=\"#eee\"/><text x=\"%.1f\" y=\"%d\" text-anchor=\"middle\">%d</text>\n", x, x, T, T + ph, x, T + ph + 16, v }
  printf "<text x=\"%d\" y=\"%d\" text-anchor=\"middle\">%s (log scale)</text>\n", L + pw / 2, H - 18, XLABEL
  printf "<text transform=\"translate(18,%d) rotate(-90)\" text-anchor=\"middle\">ours / best external (log scale)</text>\n", T + ph / 2 }
$1 == OP { best = ($5 < $6) ? $5 : $6; r = $4 / best
  x = L + pw * (log($3) - x0) / (x1 - x0); y = T + ph * (1 - (log(r) - y0) / (y1 - y0))
  pts[$2] = pts[$2] sprintf("%.1f,%.1f ", x, y) }
END { ly = 58; k = split("p61 p127 p255 p511", ks, " ")
  for (i = 1; i <= k; i++) { p = ks[i]
    printf "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2\" points=\"%s\"/>\n", col[p], pts[p]
    q = split(pts[p], pp, " ")
    for (j = 1; j <= q; j++) { if (pp[j] == "") continue; split(pp[j], c, ","); printf "<circle cx=\"%s\" cy=\"%s\" r=\"3\" fill=\"%s\"/>\n", c[1], c[2], col[p] }
    printf "<line x1=\"620\" x2=\"645\" y1=\"%d\" y2=\"%d\" stroke=\"%s\" stroke-width=\"2\"/><text x=\"650\" y=\"%d\" font-size=\"11\">%s</text>\n", ly, ly, col[p], ly + 4, lab[p]; ly += 20 }
  split("best of 2 rounds per tool, same inputs,|identical outputs; 1 isolated CPU,|4-vCPU Xeon VM, single-threaded|data: " (DATA == "" ? "results/p15-external/" : DATA) "", notes, "|")
  for (i = 1; i <= 4; i++) printf "<text x=\"620\" y=\"%d\" font-size=\"10\" fill=\"#555\">%s</text>\n", ly + 12 + 14 * (i - 1), notes[i]
  print "</svg>" }
