from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
                                PageBreak, KeepTogether, Preformatted)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

FD = "/usr/share/fonts/truetype/dejavu/"
pdfmetrics.registerFont(TTFont("DJ", FD + "DejaVuSerif.ttf"))
pdfmetrics.registerFont(TTFont("DJ-B", FD + "DejaVuSerif-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DJ-I", FD + "DejaVuSans.ttf"))       # no serif italic shipped; sans as 'italic'
pdfmetrics.registerFont(TTFont("DJ-BI", FD + "DejaVuSans-Bold.ttf"))
pdfmetrics.registerFont(TTFont("DJM", FD + "DejaVuSansMono.ttf"))
addMapping("DJ", 0, 0, "DJ"); addMapping("DJ", 1, 0, "DJ-B"); addMapping("DJ", 0, 1, "DJ-I"); addMapping("DJ", 1, 1, "DJ-BI")

body = ParagraphStyle("body", fontName="DJ", fontSize=9.6, leading=13.2, alignment=TA_JUSTIFY, spaceAfter=5)
small = ParagraphStyle("small", parent=body, fontSize=8.2, leading=10.5)
h1 = ParagraphStyle("h1", fontName="DJ-B", fontSize=13.5, leading=17, spaceBefore=12, spaceAfter=5)
h2 = ParagraphStyle("h2", fontName="DJ-B", fontSize=10.8, leading=14, spaceBefore=8, spaceAfter=3)
title = ParagraphStyle("title", fontName="DJ-B", fontSize=17, leading=21, alignment=TA_CENTER, spaceAfter=6)
sub = ParagraphStyle("sub", fontName="DJ", fontSize=9.5, leading=12, alignment=TA_CENTER, spaceAfter=3)
box = ParagraphStyle("box", parent=body, fontSize=9, leading=12, backColor=colors.Color(0.97, 0.95, 0.88),
                     borderColor=colors.Color(0.6, 0.5, 0.2), borderWidth=0.8, borderPadding=6, spaceBefore=6, spaceAfter=10)
mono = ParagraphStyle("mono", fontName="DJM", fontSize=8, leading=10, leftIndent=10, spaceAfter=5)
bul = ParagraphStyle("bul", parent=body, leftIndent=12, bulletIndent=2, spaceAfter=3)
cell = ParagraphStyle("cell", fontName="DJ", fontSize=8, leading=10)
cellb = ParagraphStyle("cellb", fontName="DJ-B", fontSize=8, leading=10)

def P(t, s=body): return Paragraph(t, s)
def B(t): return Paragraph(t, bul, bulletText="•")
def tbl(rows, widths, head=True):
    data = [[Paragraph(c, cellb if (head and i == 0) else cell) for c in r] for i, r in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1 if head else 0)
    st = [("GRID", (0, 0), (-1, -1), 0.4, colors.grey), ("VALIGN", (0, 0), (-1, -1), "TOP"),
          ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4)]
    if head:
        st.append(("BACKGROUND", (0, 0), (-1, 0), colors.Color(0.9, 0.9, 0.9)))
    t.setStyle(TableStyle(st)); t.spaceAfter = 8
    return t

S = []
S += [P("Eigen-Séta: Public-Key Encryption from a Secret Endomorphism of the Public Curve", title),
      P("Design note and toy-prototype report", sub),
      P("Research branch <font name='DJM' size='8.5'>research/endring-pke-strategy-20261010</font> of the <i>crypto</i> repository · 11 October 2026", sub),
      Spacer(1, 6),
      P("<b>Status.</b> This document describes a <i>proposed</i> scheme, a SageMath prototype that runs its complete "
        "KeyGen / Encaps / Decaps path at toy parameters, and the security analysis available today. "
        "Security is <b>not established</b>: the scheme rests on two assumptions that have no reduction to a studied problem, "
        "the prototype exercises only special-key variants at 32–66-bit primes, and all parameter figures are derived "
        "estimates, not measurements. Nothing here should be used to protect data.", box)]

S += [P("Abstract", h1),
      P("We propose a public-key encryption scheme whose decryption key is a secret endomorphism θ of the public supersingular "
        "curve E<sub>A</sub>, i.e. a piece of End(E<sub>A</sub>). The public key is E<sub>A</sub> together with an eigenbasis (P, Q) of θ on the "
        "2<super>e</super>-torsion; the ciphertext is the codomain E<sub>C</sub> of a secret 3<super>f</super>-isogeny φ from E<sub>A</sub> plus the "
        "torsion images φ(P), φ(Q) masked by a <i>sender-random diagonal</i> matrix. Decryption uses the fact, from Castryck–Vercauteren's "
        "generalised lollipop attack, that a mask cancels in the endomorphism φθφ̂ of E<sub>C</sub> exactly when it commutes with θ on the "
        "torsion: in an eigenbasis every diagonal mask does. The receiver therefore learns φθφ̂ on E<sub>C</sub>[2<super>e</super>] for free, "
        "turns it into an explicit isogeny with Kani's lemma, and reads off ker φ̂. Compared with SILBE (SAC 2024), the only published "
        "encryption scheme built on lollipop decoding, the move from scalar to diagonal masks removes the requirement that the torsion "
        "order have at least 2λ prime factors, and the derived prime size drops from 13 013 bits to about 772 bits at λ = 128 "
        "(assuming the two-halves Kani evaluation). We report a prototype, the corrections it forced on the design, a first "
        "cryptanalytic experiment, and the open problems.")]

S += [P("1. Problem and prior art", h1),
      P("<b>Endomorphism-ring hardness.</b> For a supersingular curve E/F<sub>p²</sub>, computing End(E) is equivalent, under probabilistic "
        "polynomial-time reductions and now unconditionally, to the isogeny-path, maximal-order, one-endomorphism and arbitrary-degree "
        "isogeny problems (Wesolowski 2021; Page–Wesolowski 2024; Herlédan Le Merdy–Wesolowski 2025). Worst-case hardness implies "
        "average-case hardness. Best attacks: classical Õ(p<super>1/2</super>) with low memory, heuristic p<super>1/3+o(1)</super> with exponential "
        "memory (Wesolowski 2026), quantum Õ(p<super>1/4</super>). Any <i>structure</i> revealed about End(E) collapses this: one endomorphism "
        "gives |disc|<super>1/4</super> classical and subexponential quantum; a full-rank suborder or an evaluable isogeny from a known-End curve gives "
        "quantum polynomial time (Chen et al. 2023; Chen–Petit 2025)."),
      P("<b>The structural obstruction.</b> Given End(E<sub>A</sub>) alone, computing an isogeny to a fresh curve E<sub>C</sub> is equivalent to "
        "computing End(E<sub>C</sub>). A receiver holding End(E<sub>A</sub>) therefore cannot decrypt “here is a random curve isogenous to yours”. "
        "A ciphertext must carry a <i>hint</i>, and the design question is which hint is useless to everyone, decodable with End(E<sub>A</sub>), "
        "and cheap. The seven-open-problems paper of Castryck, De Feo, Galbraith, Kutas, Reijnders and Wesolowski (ePrint 2026/1431) lists "
        "encryption from a small variant of EndRing as open Problem 7."),
      P("<b>Two known mechanisms.</b> Commutative group actions (CSIDH, SCALLOP, PEGASIS, CORAL) need a public orientation so the sender can act, "
        "which exposes a subexponential quantum attack. Torsion hints (the SIDH family: FESTA, QFESTA, POKÉ, PIKE, LIT-SiGamal, Séta, SILBE) "
        "reveal torsion images hidden by a subgroup Γ of GL<sub>2</sub>(Z/N). Unmasked or unipotent hints are broken for everyone once "
        "N² > deg φ (Castryck–Decru, Maino et al., Robert 2022). Scalar masks (M-SIDH) and diagonal masks (FESTA's CIST problem) only have "
        "exponential attacks on generic curves (De Feo–Fouotsa–Panny 2024)."),
      tbl([["Scheme", "Trapdoor", "Hint in ciphertext", "Why it is where it is"],
           ["Séta (AC 2021)", "curve with a secret endomorphism; Petit's attack as decryptor", "exact torsion images", "broken 2022: Robert's attack needs no endomorphism"],
           ["M-SIDH (EC 2023)", "key exchange", "scalar-masked images, N with ≥ 2λ prime factors", "pairing leaks λ²; needs E<sub>0</sub> of unknown End; 5911-bit p"],
           ["Castryck–Vercauteren (AC 2023)", "attack", "lollipop φ′ωφ̂; mask cancels iff ω commutes with Γ", "diagonal masks: only from an eigenbasis; treated as a backdoor"],
           ["SILBE (SAC 2024)", "End(E<sub>A</sub>) via a path from F<sub>p</sub>-rational E<sub>0</sub>", "scalar-masked images (mask = message)", "OW-PCA ⇔ M-SIDH problem; p = 13 013 bits; impractical"],
           ["FESTA / QFESTA / POKÉ", "secret isogeny from known-End E<sub>0</sub> + masks", "diagonal or scalar masks", "key recovery is a torsion-hint problem, not EndRing"],
           ["Eigen-Séta (this note)", "secret θ ∈ End(E<sub>A</sub>) with public eigenbasis", "sender-random diagonal masks in that basis", "lollipop decodable by the owner; CV's backdoor used legitimately"]],
          [3.0*cm, 4.3*cm, 4.2*cm, 5.5*cm])]

S += [P("2. The decoding lemma", h1),
      P("Let φ: E<sub>A</sub> → E<sub>C</sub> have degree d, let (P, Q) be a basis of E<sub>A</sub>[N], and let the ciphertext reveal "
        "Λ = φ∘D on E<sub>A</sub>[N], where D ∈ Γ ≤ GL<sub>2</sub>(Z/N) is the sender's mask written in the basis (P, Q). For θ ∈ End(E<sub>A</sub>) "
        "with matrix M on E<sub>A</sub>[N], the receiver can compute Λ∘M∘(d·Λ<super>−1</super>) = φ D M D<super>−1</super> φ̂ on E<sub>C</sub>[N]. "
        "This equals the lollipop ψ<sub>θ</sub> := φθφ̂ exactly iff DM = MD (Castryck–Vercauteren Lemma 3 with σ<sub>0</sub> = id). Hence:"),
      B("Γ = scalars: every θ works (M-SIDH, SILBE), but the Weil pairing leaks the square of the mask, leaving 2<super>ω(N)</super> candidates, so N needs at least 2λ prime factors."),
      B("Γ = diagonal in (P, Q): works iff M is diagonal, i.e. (P, Q) is an eigenbasis of θ. The centraliser of a non-scalar θ mod N is exactly the diagonal torus, so diagonal masks are the <b>largest</b> mask set compatible with a single secret endomorphism."),
      B("Once ψ<sub>θ</sub> is known on E<sub>C</sub>[N] with N² > deg ψ<sub>θ</sub> = d²·n(θ), Kani's lemma gives ψ<sub>θ</sub> as an explicit isogeny and ker φ̂ = ker ψ<sub>θ</sub> ∩ E<sub>C</sub>[d] in the generic case."),
      P("<b>Parity (found while implementing).</b> Eigenlines of full order 2<super>e</super> exist iff 2 splits in Z[θ], i.e. disc Z[θ] ≡ 1 mod 8, "
        "i.e. trace t odd and norm w even. Then one eigenvalue is even, ψ<sub>θ</sub> + [c] has even degree for every c, and the Kani kernel degenerates "
        "on the 2-torsion. The fix is to decode with the trace-zero element θ′ := 2θ − t, of norm Δ′ := 4w − t² ≡ 7 mod 8, which has the same "
        "eigenlines and eigenvalues ±s<sub>0</sub> with s<sub>0</sub>² ≡ −Δ′ mod 2<super>e+2</super>. Then ψ′ := φθ′φ̂ has odd eigenvalues ±d·s<sub>0</sub>.")]

S += [P("3. The scheme", h1),
      P("<b>Parameters.</b> p = 2<super>e</super>·3<super>f</super>·c − 1, N = 2<super>e</super>, d = 3<super>f</super> ≈ 2<super>2λ</super>, e fixed by §5. All arithmetic is over F<sub>p²</sub>."),
      P("KeyGen", h2),
      B("Walk from E<sub>0</sub> (j = 1728) by a secret isogeny of degree > 2<super>4λ</super> to a statistically uniform E<sub>A</sub>; keep the connecting ideal, so O<sub>A</sub> = End(E<sub>A</sub>) is known (Deuring tooling: “Deuring for the people”, SQIsign2D ideal-to-isogeny)."),
      B("Enumerate short non-scalar θ ∈ O<sub>A</sub> (LLL on the trace-zero lattice; shortest norm ≈ p<super>2/3</super> for a random maximal order). Pick θ with (i) Δ′ ≡ 7 mod 8 and (ii) N² − d²Δ′ = c² + c′² a sum of two squares (c even, c′ odd; e.g. a prime ≡ 1 mod 4, Cornacchia)."),
      B("Compute the eigenbasis (P, Q) of θ on E<sub>A</sub>[2<super>e</super>] by evaluating θ on a basis and solving the 2×2 eigenproblem mod 2<super>e</super>. Also fix a canonical basis (P<sub>d</sub>, Q<sub>d</sub>) of E<sub>A</sub>[3<super>f</super>]."),
      B("pk = (E<sub>A</sub>, P, Q). sk = (s<sub>0</sub>, c, c′). Note that θ itself is not needed at decryption time."),
      P("Encaps", h2),
      B("Sample s ∈ Z/3<super>f</super> and odd nonces a, b ∈ Z/2<super>e</super> from a seed. K := P<sub>d</sub> + [s]Q<sub>d</sub>; φ: E<sub>A</sub> → E<sub>C</sub> := E<sub>A</sub>/⟨K⟩ (3-isogeny chain)."),
      B("ct = (E<sub>C</sub>, S := [a]φ(P), T := [b]φ(Q)); key = H(seed, ct)."),
      P("Decaps", h2),
      B("ψ′ acts on the ciphertext points by known scalars: ψ′(S) = [d·s<sub>0</sub>]S, ψ′(T) = [−d·s<sub>0</sub>]T, and ψ̂′ = −ψ′. (Derivation: ψ′(aφP) = aφθ′(φ̂φP) = a·d·φ(θ′P) = a·d·s<sub>0</sub>·φ(P); the nonces never enter.)"),
      B("Kani in dimension 2: with α := ψ′ + [c] and c² + c′² = N² − d²Δ′, the matrix F = [[α, c′], [−c′, α̂]] is an (N², N²)-isogeny of E<sub>C</sub>² with kernel {(c′x, −αx) : x ∈ E<sub>C</sub>[N²]}, whose N-torsion part ⟨(c′S, −(ds<sub>0</sub>+c)S), (c′T, −(c−ds<sub>0</sub>)T)⟩ is known. Compute F as two N-halves from the actions of α and α̂; then α(X) is the first component of F(X, 0) and c′X the second, which also fixes the component isomorphisms of the split codomain."),
      B("ker φ̂ = ker(α − c) ∩ E<sub>C</sub>[3<super>f</super>] (cyclic of order d generically; otherwise try the few cyclic order-d subgroups). Recover φ, K and s; re-derive a, b from the seed, re-encrypt and compare (Fujisaki–Okamoto with implicit rejection)."),
      P("<b>Correctness.</b> ψ′ = φθ′φ̂ is an endomorphism of E<sub>C</sub> of degree d²Δ′ and trace 0, so deg α = d²Δ′ + c² and deg α + c′² = N². "
        "Maximal isotropy of the kernel follows from e(c′S, c′T)·e(αS, αT) = e(S, T)<super>c′² + deg α</super> = 1. Because θ′ is diagonal with odd eigenvalues it acts as the "
        "identity on E<sub>C</sub>[2]; the kernel's 2-torsion is the diagonal of E<sub>C</sub>[2]², so the (2,2)-chain necessarily starts with product-to-product steps "
        "before gluing into a Jacobian. Generic ker φ ⊄ θ′-stable subgroups gives ker ψ′ ∩ E<sub>C</sub>[d] = ker φ̂.")]

S += [PageBreak(), P("4. Prototype", h1),
      P("A SageMath implementation of the full path lives in <font name='DJM' size='8.5'>research/eigen_seta_prototype_20261011/</font>. Because the "
        "FESTA and Castryck–Decru chain code assumes one gluing, Jacobian steps, and one split, a generic (2<super>e</super>, 2<super>e</super>)-chain evaluator "
        "(<font name='DJM' size='8.5'>kani_chain.py</font>) was written that handles product→product, gluing, Richelot and splitting steps in any order, "
        "reusing the Richelot/Mumford formulas of Pope–Oudompheng–Panny (MIT)."),
      P("<b>Two special-key variants are implemented.</b> <i>cm</i>: E<sub>A</sub> is the curve with CM by Z[(1+√−7)/2] (j = −3375, supersingular for p ≡ 3, 5, 6 mod 7); "
        "θ = ω of degree 2 (t = ±1, w = 2, Δ′ = 7). Its End is public, so it is insecure as a key but exercises the whole mechanism. "
        "<i>walk</i>: E<sub>A</sub> = ε(E<sub>cm</sub>) for a secret 5<super>k</super>-isogeny ε and θ = εωε̂ (Δ′ = 7·25<super>k</super>): Castryck–Vercauteren's §5.3 backdoor "
        "used by its owner. Both run the single-chain regime (N > deg ψ′). The random-curve variant of §3 is <b>not</b> implemented: it needs the two-halves evaluation (§5)."),
      P("Results (toy sizes, correctness only)", h2),
      tbl([["variant", "p (bits)", "chain length e", "sender degree", "trials", "correct", "decaps time (Sage)"],
           ["cm", "41", "26", "3<super>5</super>", "25", "25", "≈ 0.2 s"],
           ["walk, 5<super>2</super>", "45", "29", "3<super>4</super>", "25", "25", "≈ 0.2 s"],
           ["walk, 5<super>2</super>", "55", "36", "3<super>5</super>", "25", "25", "≈ 0.3 s"],
           ["walk, 5<super>3</super>", "66", "44", "3<super>6</super>", "10", "10", "≈ 0.8 s"]],
          [2.6*cm, 1.8*cm, 2.6*cm, 2.6*cm, 1.6*cm, 1.8*cm, 3.6*cm]),
      P("In every trial the independent check of the lemma (direct evaluation of φθ′φ̂ on the ciphertext points) matched the predicted eigenvalues ±d·s<sub>0</sub>, "
        "and the chain began with two or three product-type steps before gluing, as predicted. Two corner cases were met and handled: a non-cyclic "
        "ker ψ′ ∩ E<sub>C</sub>[3<super>f</super>] (resolved by testing the cyclic order-d candidates against the codomain E<sub>A</sub>) and a degenerate Mumford input in the "
        "splitting formulas for some payload points (resolved by retrying with a fresh basis of E<sub>C</sub>[3<super>f</super>]; about 1 in 25 trials; cause not isolated)."),
      P("Experiment: do wrong eigenvalue guesses fail early?", h2),
      P("An attacker who tries to recover θ from the public eigenbasis must guess Δ′ (hence s<sub>0</sub>, c, c′). Every pairing-consistent guess gives a maximal isotropic "
        "kernel, so the chain runs; the question is whether a wrong guess fails only at the final split or earlier (which would allow a digit-by-digit attack, as in "
        "Castryck–Decru). At e = 36 (walk) and e = 26 (cm), with 43 guesses each including guesses agreeing with the true Δ′ on 3–13 low bits:"),
      tbl([["outcome", "walk, e = 36", "cm, e = 26"],
           ["wrong guess ran the full chain and ended on a Jacobian (no split)", "43 / 43", "42 / 43"],
           ["wrong guess split to E<sub>C</sub> × E<sub>C</sub> (false positive)", "0", "1 (random guess)"],
           ["true value split to E<sub>C</sub> × E<sub>C</sub>", "yes", "yes"]],
          [10.5*cm, 3.2*cm, 3.2*cm]),
      P("No early-abort signal appeared. The single false positive is explained by the CM toy curve: E<sub>C</sub> is 3<super>f</super>-isogenous to a curve of discriminant −7, "
        "so other elements of the same quadratic order share the eigenbasis and some wrong Δ′ are genuine. This is an artefact of the special key and, at the same time, "
        "a warning for special-key variants in general.")]

S += [P("5. Parameters and the two-halves requirement", h1),
      P("Constraints: meet-in-the-middle on the 3<super>f</super>-isogeny forces d ≥ 2<super>2λ</super>; decoding needs N > d√Δ′; for a random curve the shortest "
        "odd-trace endomorphism has norm w ≈ p<super>2/3</super>; and p ≈ N·d. Solving gives N ≈ d², p ≈ d³. Derived figures at λ = 128:"),
      tbl([["quantity", "value", "comparison"],
           ["sender degree d = 3<super>f</super>", "3<super>162</super> ≈ 2<super>257</super>", "FESTA uses 257–260-bit degrees"],
           ["eigenbasis torsion N = 2<super>e</super>", "e ≈ 515", "N > d·p<super>1/3</super>"],
           ["prime p", "≈ 2<super>772</super>", "FESTA 1292, POKÉ 431, SILBE 13 013 bits"],
           ["public key", "≈ 325 bytes (curve + compressed eigenbasis)", "FESTA 561, POKÉ 324"],
           ["ciphertext", "≈ 325 bytes (curve + compressed masked points)", "FESTA 1122, POKÉ 648"],
           ["decryption", "one (2<super>515</super>, 2<super>515</super>)-chain in dimension 2 as two halves + 3<super>162</super> kernel recovery", "FESTA-class cost"]],
          [4.2*cm, 6.4*cm, 6.4*cm]),
      P("<b>Why the two halves are essential.</b> The table assumes the (N², N²)-isogeny is computed from N-torsion as two N-halves (N² > deg). With a single chain "
        "(N > deg, the regime the prototype implements) the condition becomes N > d²Δ′ ≈ 4d²p<super>2/3</super> with N ≤ p, which forces p > d<super>9</super> ≈ 2<super>18λ</super>; "
        "and for curves close to E<sub>0</sub> every odd-trace endomorphism has norm at least p/(4D<sub>A</sub>), which the prototype confirmed numerically, so a single chain never "
        "fits a random key. In Richelot/Mumford form the second half needs the dual chain plus an explicit isomorphism of the middle Jacobians; in theta coordinates "
        "(SQIsignHD, SQIsign2D) duals are explicit. This is the next implementation milestone."),
      P("<b>Hidden-CM-curve variant at full size.</b> Finding the secret walk ε costs √D<sub>A</sub> and finding εωε̂ costs about √7·D<sub>A</sub>, so D<sub>A</sub> ≈ 2<super>2λ</super>, "
        "Δ′ ≈ 2<super>4λ</super>, N ≈ 2<super>4λ</super>, p ≈ 2<super>6λ</super>: the same prime size as the random-curve variant, a far cheaper key generation, "
        "but a different and probably weaker key-recovery assumption (a hidden small-discriminant orientation plus a walk, rather than EndRing).")]

S += [P("6. Security analysis", h1),
      P("Two assumptions are new and have no reduction to a studied problem."),
      B("<b>Eigen-CIST</b> (message recovery): given (E<sub>A</sub>, P, Q, E<sub>C</sub>, [a]φP, [b]φQ) with random odd a, b, find φ. The Weil pairing leaks ab·d mod 2<super>e</super>, "
        "leaving ≈ 2<super>e−1</super> masks; an attacker needs the correct mask on a block 2<super>e′</super> with 2<super>2e′</super> > d before Robert's attack applies, cost ≈ √d ≈ 2<super>λ</super>. "
        "This is FESTA's CIST problem family, unbroken since 2023 on generic bases, except that here the basis is a secret eigenbasis."),
      B("<b>EndRing with an eigenbasis hint</b> (key recovery): given (E<sub>A</sub>, P, Q), find θ. Guessing the eigenvalues modulo a block 2<super>e′</super> with 2<super>2e′</super> > Δ′ and running "
        "Kani costs ≈ Δ′ ≈ 2<super>2λ</super>; cycle search for an endomorphism of norm Δ′ costs ≈ 2<super>2λ</super>; generic EndRing costs p<super>1/2</super> ≈ 2<super>386</super> "
        "(p<super>1/3</super> ≈ 2<super>257</super> heuristically). Once θ is known, End(E<sub>A</sub>) follows at |disc|<super>1/4</super> classical or subexponential quantum cost "
        "(Herlédan Le Merdy–Wesolowski). The eigenlines ⟨P⟩, ⟨Q⟩ also expose two horizontal 2-walks of length e in the Z[θ]-orientation, a negligible fraction of the class group; "
        "self-pairing attacks (Castryck et al. 2023) need m² | Δ′ with m | N and are excluded by Δ′ odd."),
      tbl([["attack", "best known cost at the §5 parameters", "status"],
           ["mask guessing + Robert (ciphertext)", "≈ 2<super>128</super> classical; ≈ 2<super>64</super> with Grover", "same as FESTA and every masked-torsion scheme"],
           ["meet-in-the-middle on the 3<super>162</super>-isogeny", "≈ 2<super>128</super>", "generic"],
           ["eigenvalue guessing + Kani (key)", "≈ 2<super>256</super>", "no early-abort signal at toy size (§4)"],
           ["cycle search for θ", "≈ 2<super>256</super>", "generic"],
           ["generic EndRing", "2<super>386</super> / 2<super>257</super> / quantum 2<super>193</super>", "not the binding constraint"],
           ["hidden-CM variant: find the walk ε", "≈ √D<sub>A</sub>", "sets D<sub>A</sub> ≈ 2<super>2λ</super>; weaker assumption than EndRing"]],
          [5.4*cm, 5.6*cm, 6.0*cm]),
      P("<b>Level achieved.</b> The prototype provides no security: its primes are 32–66 bits, the cm key has public End, and the walk keys are found by exhaustive search in seconds. "
        "The design at the §5 parameters has a <i>conjectured</i> classical level of 128 bits and, like every masked-torsion scheme, about 64 bits against a quantum mask-guessing "
        "attack unless d is raised to 2<super>4λ</super>. The figure is a conjecture because the two assumptions above are unproven; the only evidence is the toy experiment of §4."),
      P("<b>Other requirements.</b> Decryption uses the secret s<sub>0</sub>, so a chosen-ciphertext attacker with a decryption oracle learns about it unless the Fujisaki–Okamoto "
        "transform with re-encryption and implicit rejection is in place (compare the Moriya–Onuki adaptive attacks on FESTA). Masks must be fresh per ciphertext and derived from "
        "the seed; partial mask leakage (Das–Invernizzi–Kutas–Meers 2025) and fault attacks on torsion masking (Gilchrist–Lai–Meyer 2026) apply as to FESTA. If θ ever leaks "
        "through a side channel, security drops to the subexponential one-endomorphism bound: constant-time Kani code is mandatory from the start.")]

S += [P("7. What is new, and what is not", h1),
      B("vs Séta: masks make the hint useless to Robert's generic attack; torsion requirement N > d√Δ′ instead of N > D²; no trapdoor-curve construction."),
      B("vs SILBE: diagonal instead of scalar masks removes the ω(N) ≥ 2λ constraint (13 013 → ≈ 772 bits); no F<sub>p</sub>-rational E<sub>0</sub>, no Frobenius; a shortest endomorphism replaces the degree-√p connecting isogeny. Lost: SILBE's public key is any curve with no extra data; ours leaks an eigenbasis."),
      B("vs FESTA/QFESTA/POKÉ: the public key carries no torsion images of a receiver isogeny from a known-End curve, so key recovery is an EndRing-type problem rather than CIST or C-POKE; E<sub>0</sub> is not used at encryption time."),
      B("vs Castryck–Vercauteren §5.3: the “backdoored basis” is the honest receiver's secret, so the basis-hashing countermeasure does not apply."),
      P("The scheme does <b>not</b> solve open Problem 7 of ePrint 2026/1431: its ciphertext assumption is still a torsion-hint assumption. It would be the first practical-size "
        "encryption scheme whose key-recovery problem is an EndRing variant with a leakage term rather than a torsion-hint problem, if the assumptions hold.")]

S += [P("8. Open problems and next steps", h1),
      B("Implement the two-halves Kani evaluation (dual Richelot chain plus middle-Jacobian isomorphism, or a theta-model implementation) and the random-curve key generation (LLL on End(E<sub>A</sub>) with explicit evaluation of a short endomorphism)."),
      B("Repeat the wrong-guess experiment on random-curve keys, where the CM false positive cannot occur, and at larger sizes."),
      B("Study eigen-CIST against the level-structure reductions of De Feo–Fouotsa–Panny and the generalised class-group actions of Arpin et al., with the eigenbasis as the level structure."),
      B("Parameter script and a C implementation; constant-time chain code; Fujisaki–Okamoto wrapper."),
      B("Formalise Eigen-EndRing and eigen-CIST; prove OW-CPA ⇐ eigen-CIST and key recovery ⇐ Eigen-EndRing.")]

S += [P("References (ePrint identifiers)", h1),
      P("Wesolowski, FOCS 2021 (2021/919). Page–Wesolowski, EC 2024 (2023/1399). Herlédan Le Merdy–Wesolowski, TCC 2025 (2025/271); CiC 2025 (2023/1448). "
        "Wesolowski 2026/1486; Mamah 2026/1821. Castryck–De Feo–Galbraith–Kutas–Reijnders–Wesolowski, “The isogeny problems”, 2026/1431. "
        "Castryck–Vercauteren, AC 2023 (2023/1433). Fouotsa–Moriya–Petit, EC 2023 (2023/013). De Feo et al., Séta, AC 2021 (2019/1291). "
        "Duparc–Fouotsa–Vaudenay, SILBE, SAC 2024 (2024/400). Basso–Maino–Pope, FESTA, AC 2023 (2023/660). Nakagawa–Onuki, QFESTA, Crypto 2024 (2023/1468). "
        "Basso–Maino, POKÉ, EC 2025 (2024/624). Basso et al., UPKE from FESTA, AC 2026 (2026/1014). De Feo–Fouotsa–Panny, EC 2024 (2024/459). "
        "Chen–Imran–Ivanyos–Kutas–Leroux–Petit, AC 2023 (2023/779); Chen–Petit, EC 2025. Castryck–Houben–Merz–Mula–van Buuren–Vercauteren, Crypto 2023 (2023/549). "
        "Robert, EC 2023 (2022/1038). Dartois–Leroux–Robert–Wesolowski, SQIsignHD (2023/436). Basso et al., SQIsign2D-West (2024/760). "
        "Das–Invernizzi–Kutas–Meers, PKC 2026 (2025/2075). Gilchrist–Lai–Meyer, AC 2026 (2026/1944). Pope–Oudompheng–Panny, Castryck-Decru-SageMath (GitHub, MIT).", small)]

def footer(canvas, doc):
    canvas.saveState(); canvas.setFont("DJ", 7.5); canvas.setFillColor(colors.grey)
    canvas.drawString(2*cm, 1.2*cm, "Eigen-Séta design note · toy prototype · security not established")
    canvas.drawRightString(A4[0] - 2*cm, 1.2*cm, f"page {doc.page}")
    canvas.restoreState()

out = "/home/user/crypto/research/eigen_seta_prototype_20261011/eigen_seta_design_note.pdf"
doc = SimpleDocTemplate(out, pagesize=A4, leftMargin=2*cm, rightMargin=2*cm, topMargin=1.8*cm, bottomMargin=1.8*cm,
                        title="Eigen-Séta: Public-Key Encryption from a Secret Endomorphism of the Public Curve",
                        author="crypto repository research branch", subject="Isogeny-based PKE design note and prototype report")
doc.build(S, onFirstPage=footer, onLaterPages=footer)
print("wrote", out)
