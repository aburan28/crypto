import { readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const rows = readFileSync(join(here, "cm-degree-frontier.csv"), "utf8")
  .trim()
  .split("\n")
  .slice(1)
  .map((line) => {
    const [id, label, bits, status, source] = line.split(",");
    return { id, label, bits: Number(bits), status, source };
  });

const width = 1060;
const height = 470;
const left = 340;
const right = 48;
const top = 84;
const rowHeight = 50;
const axisWidth = width - left - right;
const maximum = 200;
const x = (value) => left + (value / maximum) * axisWidth;
const escapeXml = (value) => value
  .replaceAll("&", "&amp;")
  .replaceAll("<", "&lt;")
  .replaceAll(">", "&gt;")
  .replaceAll('"', "&quot;");
const color = (row) => {
  if (row.id.startsWith("ramified")) return "#7e22ce";
  if (row.id === "exact-search-boundary") return "#2563eb";
  if (row.id === "non-scalar-lower-bound") return "#b91c1c";
  return "#64748b";
};

const tickValues = [0, 48, 96, 144, 192, 200];
const ticks = tickValues.map((value) => `
  <line x1="${x(value).toFixed(2)}" y1="${top - 12}" x2="${x(value).toFixed(2)}" y2="${height - 58}" class="grid"/>
  <text x="${x(value).toFixed(2)}" y="${height - 37}" text-anchor="middle" class="tick">${value}</text>`).join("");
const bars = rows.map((row, index) => {
  const y = top + index * rowHeight;
  const barWidth = Math.max(2, x(row.bits) - left);
  return `
  <g id="${escapeXml(row.id)}">
    <text x="${left - 14}" y="${y + 18}" text-anchor="end" class="label">${escapeXml(row.label)}</text>
    <rect x="${left}" y="${y}" width="${barWidth.toFixed(2)}" height="25" rx="4" fill="${color(row)}"/>
    <text x="${Math.min(x(row.bits) + 8, width - right - 78).toFixed(2)}" y="${y + 18}" class="value">${row.bits.toFixed(3)} bits</text>
    <text x="${left - 14}" y="${y + 35}" text-anchor="end" class="status">${escapeXml(row.status)} · ${escapeXml(row.source)}</text>
  </g>`;
}).join("");

const svg = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="title description">
  <title id="title">P-192 CM exact search boundary and theorem control</title>
  <desc id="description">Derived and frozen logarithmic ideal norm values. These are not timing measurements.</desc>
  <style>
    text { font-family: "DejaVu Sans", sans-serif; fill: #0f172a; }
    .title { font-size: 19px; font-weight: 700; }
    .subtitle { font-size: 12px; fill: #475569; }
    .label { font-size: 12px; font-weight: 600; }
    .status { font-size: 9px; fill: #64748b; }
    .value { font-size: 11px; font-weight: 700; }
    .tick { font-size: 10px; fill: #475569; }
    .grid { stroke: #cbd5e1; stroke-width: 1; stroke-dasharray: 3 5; }
    .axis { stroke: #334155; stroke-width: 1.3; }
  </style>
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="32" y="30" class="title">CM degree frontier: exact scope versus non-scalar theorem bound</text>
  <text x="32" y="52" class="subtitle">Unit: log₂ ideal norm / geometric degree. Derived and preregistered values only; no run or timing result is plotted.</text>
${ticks}
  <line x1="${left}" y1="${height - 58}" x2="${width - right}" y2="${height - 58}" class="axis"/>
${bars}
  <text x="${left + axisWidth / 2}" y="${height - 10}" text-anchor="middle" class="subtitle">log₂(norm)</text>
</svg>
`;

writeFileSync(join(here, "cm-degree-frontier.svg"), svg);
