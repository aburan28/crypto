#!/bin/sh
# Editable source for oracle_ladder.svg and product_law.svg.
# Coordinates: x_px = 90 + (l-4)*152.5 ; y_px = 400 - 15*log2(count).
# Data (frozen artefacts, see the note's §2):
#   native CDCL+XOR conflicts, uniform rejection, l=5..8: 5781, 45808, 363293, 2877905
#     (research/notes/index-calculus/RESEARCH_SAT_SEMAEV.md, "Does it scale? No.")
#   CryptoMiniSat conflicts, uniform exhaustion on K0/F_2^131, d=4..8: 1, 1, 1, 839, 312231
#     (research/nagao_relations/solver_16/RESULTS.md)
#   plain WDSat decisions = 2^{3d}/3! to 0.2 % (same source)
#   enumeration floor T_enum = C(2^l, 2) pairs; full-tuple line 2^{3l}/3!; MITM 2^l lookups.
set -e
cat > oracle_ladder.svg <<'SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="470" viewBox="0 0 1000 470" font-family="Helvetica, Arial, sans-serif" font-size="13">
  <rect width="1000" height="470" fill="#ffffff"/>
  <text x="500" y="24" text-anchor="middle" font-size="16" font-weight="bold">Per-target oracle work on the Weil-descended S4 system, m = 3</text>
  <text x="500" y="42" text-anchor="middle" fill="#555">solver counter per uniform target (log2) against the enumeration floor; measured points, lines derived</text>
  <!-- axes -->
  <line x1="90" y1="400" x2="700" y2="400" stroke="#222" stroke-width="1.2"/>
  <line x1="90" y1="400" x2="90" y2="60" stroke="#222" stroke-width="1.2"/>
  <!-- y grid -->
  <g stroke="#ddd" stroke-width="1">
    <line x1="90" y1="340" x2="700" y2="340"/><line x1="90" y1="280" x2="700" y2="280"/><line x1="90" y1="220" x2="700" y2="220"/><line x1="90" y1="160" x2="700" y2="160"/><line x1="90" y1="100" x2="700" y2="100"/>
  </g>
  <g fill="#333" text-anchor="end">
    <text x="82" y="404">2^0</text><text x="82" y="344">2^4</text><text x="82" y="284">2^8</text><text x="82" y="224">2^12</text><text x="82" y="164">2^16</text><text x="82" y="104">2^20</text>
  </g>
  <g fill="#333" text-anchor="middle">
    <text x="90" y="420">4</text><text x="242.5" y="420">5</text><text x="395" y="420">6</text><text x="547.5" y="420">7</text><text x="700" y="420">8</text>
    <text x="395" y="445">factor-base dimension l (subspace V, |F| ≈ 2^l), n ≈ 3l or n = 131</text>
  </g>
  <text x="30" y="230" transform="rotate(-90 30 230)" text-anchor="middle" fill="#333">conflicts or decisions per target</text>
  <!-- full-tuple line 2^{3l}/6 -->
  <polyline points="90,258.7 242.5,213.9 395,168.9 547.5,123.9 700,78.9" fill="none" stroke="#b03a2e" stroke-width="1.5" stroke-dasharray="6,4"/>
  <text x="705" y="76" fill="#b03a2e">2^(3l)/3!  all m summands assigned</text>
  <!-- enumeration floor C(2^l,2) -->
  <polyline points="90,296.4 242.5,265.8 395,235.3 547.5,205.2 700,175.2" fill="none" stroke="#1f6f8b" stroke-width="2"/>
  <text x="705" y="178" fill="#1f6f8b">T_enum = C(2^l, 2): m−1 enumerated, last root-found</text>
  <!-- MITM 2^l -->
  <polyline points="90,340 242.5,325 395,310 547.5,295 700,280" fill="none" stroke="#2e7d32" stroke-width="1.5" stroke-dasharray="2,3"/>
  <text x="705" y="283" fill="#2e7d32">2^l lookups (pair table, 2^(2l) memory)</text>
  <!-- native CDCL points -->
  <g fill="#b03a2e">
    <circle cx="242.5" cy="212.5" r="5"/><circle cx="395" cy="167.8" r="5"/><circle cx="547.5" cy="122.9" r="5"/><circle cx="700" cy="78.1" r="5"/>
  </g>
  <!-- WDSat plain points (equal to the dashed line) -->
  <g fill="none" stroke="#b03a2e" stroke-width="1.5">
    <rect x="85" y="253.7" width="10" height="10"/><rect x="237.5" y="208.9" width="10" height="10"/><rect x="390" y="163.9" width="10" height="10"/><rect x="542.5" y="118.9" width="10" height="10"/><rect x="695" y="73.9" width="10" height="10"/>
  </g>
  <!-- CryptoMiniSat on K0/F_2^131 -->
  <g fill="#6a1b9a">
    <polygon points="90,392 96,402 84,402"/><polygon points="242.5,392 248.5,402 236.5,402"/><polygon points="395,392 401,402 389,402"/><polygon points="547.5,246.4 553.5,256.4 541.5,256.4"/><polygon points="700,118.3 706,128.3 694,128.3"/>
  </g>
  <text x="300" y="392" fill="#6a1b9a" font-size="12">d ≤ 6: one conflict, the linear NO-certificate of the S4 value set (saturates at d = 7)</text>
  <!-- legend -->
  <g font-size="12">
    <circle cx="110" cy="75" r="5" fill="#b03a2e"/><text x="120" y="79">native CDCL + Gauss–Jordan XOR, uniform rejection, n ≈ 3l (RESEARCH_SAT_SEMAEV.md)</text>
    <rect x="105" y="88" width="10" height="10" fill="none" stroke="#b03a2e" stroke-width="1.5"/><text x="120" y="97">WDSat, plain mode, K0 over F_2^131 (solver_16): 2^(3d)/3! to 0.2 %</text>
    <polygon points="110,104 116,114 104,114" fill="#6a1b9a"/><text x="120" y="113">CryptoMiniSat 5.14, native XOR, K0 over F_2^131 (solver_16); d = 9, 10 censored at 300 s</text>
  </g>
  <text x="500" y="462" text-anchor="middle" fill="#555" font-size="11">Every measured point sits on or above the floor; the SAT arms sit 2^l/3 above it. Lines are derived counts, not fits. Population: 4–16 targets per cell.</text>
</svg>
SVG
cat > product_law.svg <<'SVG'
<svg xmlns="http://www.w3.org/2000/svg" width="820" height="330" viewBox="0 0 820 330" font-family="Helvetica, Arial, sans-serif" font-size="13">
  <rect width="820" height="330" fill="#ffffff"/>
  <text x="410" y="26" text-anchor="middle" font-size="16" font-weight="bold">Where a solver constant sits in the index-calculus product law (ECC2K-130, Frobenius-stable base)</text>
  <!-- three boxes -->
  <g stroke="#1f6f8b" stroke-width="1.5" fill="#eef6f9">
    <rect x="30" y="60" width="200" height="90" rx="6"/><rect x="290" y="60" width="200" height="90" rx="6"/><rect x="550" y="60" width="240" height="90" rx="6"/>
  </g>
  <g text-anchor="middle" fill="#123">
    <text x="130" y="84" font-weight="bold">relations needed</text><text x="130" y="106">|F| / n = 2^l / 131</text><text x="130" y="128" fill="#555" font-size="12">Frobenius collapse: ÷ n</text>
    <text x="390" y="84" font-weight="bold">targets per relation</text><text x="390" y="106">2^131 / C(2^l, m)</text><text x="390" y="128" fill="#555" font-size="12">yield law, measured on 12 toy cells</text>
    <text x="670" y="84" font-weight="bold">oracle per target</text><text x="670" y="106">c · C(2^l, m−1)</text><text x="670" y="128" fill="#555" font-size="12">c = solver constant vs. enumeration</text>
  </g>
  <text x="260" y="110" text-anchor="middle" font-size="22">×</text><text x="520" y="110" text-anchor="middle" font-size="22">×</text>
  <!-- result -->
  <text x="410" y="190" text-anchor="middle" font-size="15">=  c · m · 2^131 / 131   (every 2^l cancels: l, m and the base move work between boxes, never the total)</text>
  <line x1="60" y1="205" x2="760" y2="205" stroke="#999"/>
  <g font-size="13">
    <text x="60" y="230">Pollard rho with ⟨−1⟩ × ⟨π⟩ (reference):</text><text x="430" y="230" font-family="monospace">2^60.81</text>
    <text x="60" y="252">c needed to reach rho (m = 3):</text><text x="430" y="252" font-family="monospace" fill="#b03a2e">2^−64.8</text>
    <text x="60" y="274">c measured for every SAT / SMT / WDSat arm in the tree:</text><text x="430" y="274" font-family="monospace" fill="#b03a2e">2^l / m  ≥  1</text>
    <text x="60" y="296">c for pairs-and-solve, pair table, MITM:</text><text x="430" y="296" font-family="monospace">1  (2^−l with 2^(2l) memory)</text>
  </g>
  <text x="410" y="320" text-anchor="middle" fill="#555" font-size="11">Sources: RESEARCH_ECC2K130_DECOMPOSITION.md §5.1–5.3 and §6 (derived), the §2 table of this note (measured). Nothing here is a run.</text>
</svg>
SVG
echo ok
