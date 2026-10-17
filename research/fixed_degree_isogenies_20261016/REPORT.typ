// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]
#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

#let conf(
  title: none,
  authors: none,
  date: none,
  abstract: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "a4",
  lang: "en",
  region: "US",
  font: ("DejaVu Serif",),
  fontsize: 10pt,
  sectionnumbering: none,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: "1",
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)

  if title != none {
    align(center)[#block(inset: 2em)[
      #text(weight: "bold", size: 1.5em)[#title]
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[Abstract] #h(1em) #abstract
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}
#show: doc => conf(
  title: [Fixed-degree isogenies between supersingular curves, and the rest of √élu: implementation and measurements],
  date: [2026-10-16],
  cols: 1,
  doc,
)


= Fixed-degree isogenies between supersingular curves, and the rest of √élu
#strong[Question.] Are all algorithms in two sources implemented in `isogeny_algos`, and if not, implement them end to end?

+ Benčina, Kutas, Merz, Petit, Stopar and Weitkämper, #emph[Finding isogenies of fixed degree between supersingular elliptic curves] (Isogeny Club, 12 December 2023; the slides were attached to the request).
+ Bernstein, De Feo, Leroux and Smith, #emph[Faster computation of isogenies of large prime degree] (√élu, 2020; PDF attached).

#strong[Answer.] Before this change, only √élu’s core (Algorithm 2: the Weierstrass and Montgomery codomain and x-map) existed. Now every algorithm in both sources is implemented, tested against an independent computation, and run end to end, with these exceptions:

- the quantum speed-ups (Grover, Tani), for which there is only a query-count model;
- the rigorous Bauer–Joux step of the trivariate method, which is replaced by a run-time independence check.

All limits are listed in "Limits and negative results".

- #strong[Status:] implemented and measured (this run), with the scope stated per item. Every claim below is either a test in the crate or a value in the run’s records.
- #strong[Run:] `run-20261016T0959Z`, commit `06e5deb47`; records in #link("data/records.jsonl")[`data/records.jsonl`]; host in #link("data/host.txt")[`data/host.txt`].
- #strong[Host and conditions:] a 2-vCPU Xeon \@ 2.8 GHz cloud container, shared with the build. One run per setting, so the wall times carry no interval (L0).

#figure([#image("figures/fig5_pipeline_diagram.png")],
  caption: [
    The six steps as implemented
  ]
)

#emph[The six steps of the slides as implemented. Steps 1 and 6 use the Deuring correspondence on E₀ over F\_{p⁴}, which limits p to below 2³¹. Steps 2–5 are pure quaternion arithmetic and run at any size.]

== What was implemented
The code is in `isogeny_algos/src/`; the full status table is in `isogeny_algos/docs/SURVEY.md` §P1, P4 and P7 (the "V4" entries).

=== Benčina–Kutas–Merz–Petit–Stopar–Weitkämper (the slides’ steps 1–6)
#strong[Step 1. End(E₁), End(E₂).] In `fixed_degree/endring.rs`.

- #emph[Algorithm:] a bidirectional BFS on the Φ₂ j-graph from 1728 to j(E), costing Õ(√p) vertices. The path is then walked on the quaternion side: the l + 1 sub-ideals I·l + I·β (β ∈ O\_R(I)) are exactly the l-isogenies, and the one whose Deuring image has the next j is kept. The ideal is reduced to minimal norm after each step. This is the classical form of the Eisenträger–Hallgren–Lauter–Morrison–Petit reduction.
- #emph[Independent check:] the Deuring image of the result is j(E). The l + 1 sub-ideals’ images are exactly the roots of Φ\_l(j, Y), checked for l \= 2 and 3.

#strong[Step 2. Connecting ideal.] In `fixed_degree/connecting.rs`.

- #emph[Algorithm:] Kirschmer–Voight, I \= N·O₁O₂.
- #emph[Independent check:] O\_L(I) \= O₁ and O\_R(I) \= O₂. It also matches conj(I₁)·I₂, and the Deuring image of I₁I is E₂.

#strong[Step 3. Norm form.] In `connecting::norm_form`.

- #emph[Algorithm:] Q \= 2·Nrd/N(I) on an LLL-reduced basis, so the Gram matrix is integral.
- #emph[Independent check:] Q(x)·N(I) \= 2·Nrd(Σ xᵢσᵢ) for random x.

#strong[Step 4a. Cornacchia.] In `fixed_degree/cornacchia.rs` and `factor.rs`.

- #emph[Algorithm:]
  - guess (x₃, x₄);
  - complete the square to X² + |D|Y² \= N;
  - factor N, with Pollard rho;
  - take every √−D mod N by CRT;
  - abandon the guess when ω(N) \> 11·log log N.
- #emph[Independent check:] every solution set equals brute force. This covers x² + dy² \= m for all m ≤ 2500 and 12 values of d, and 300 random binary quadratic equations.

#strong[Step 4b. Coron’s bivariate method.] In `fixed_degree/coppersmith.rs`, with an exact integer LLL in `lattice.rs` and multivariate integer polynomials and Bareiss resultants in `mpoly.rs`.

- #emph[Algorithm:] a box lattice; a vector is accepted only with the certificate ‖h(xX)‖₁ \< n; then Res\_y, and the integer roots are found exactly.
- #emph[Independent check:] planted roots are recovered up to Coron’s bound (reached at k \= 2).

#strong[Step 4c. Trivariate method (Bauer–Joux setting).] In `fixed_degree/coppersmith.rs`. This one is #strong[partial].

- #emph[Algorithm:] a Coron-style trivariate lattice. Instead of Bauer–Joux’s rigorous Gröbner-basis step, the algebraic independence of the two polynomials is checked at run time by resultants modulo primes.
- #emph[Independent check:] planted roots, but only for |xᵢ| ≲ 4 (30-bit Gram entries) to ≲ 12 (40-bit).

#strong[Step 5. The ideal of norm d.] In `pipeline::solve_degree`.

- #emph[Algorithm:] J \= I·conj(α)/N(I).
- #emph[Independent check:] N(J) \= d, O\_L(J) \= O₁, and O\_R(J) \= β⁻¹O₂β exactly.

#strong[Step 6a. Ideal → isogeny, smooth d.] In `fixed_degree/represent.rs`.

- #emph[Algorithm:] the kernel E₀\[K + d·O₀\] is pushed through E₀ → E₁′ and turned into a Vélu chain.
- #emph[Independent check:] the chain from the caller’s own E₁ lands on j(E₂).

#strong[Step 6b. Ideal → isogeny, any odd d.] Same file.

- #emph[Algorithm:] compute χ(S₁), χ(S₂) on E₁\[2^(a+2)\], the Robert-style representation, through φ\_K \= \[N(K′)⁻¹\]·φ\_K′ ∘ θ.
- #emph[Certificate (Kani):] take an auxiliary isogeny of degree 2ᵃ − d, from a second random ideal and evaluated the same way. The (2ᵃ, 2ᵃ)-isogeny of X × E₂ with kernel {(g(S), χ(S))} must split, and a twisted control kernel must not.
- #emph[Independent check:] the Kani split, the twisted control, e(χS₁, χS₂) \= e(S₁, S₂)^d, and the Deuring j.

#strong[Step 7. Hybrid.] In `fixed_degree/hybrid.rs`.

- #emph[Algorithm:] guess an l^{e₁} prefix as a chain of sub-ideals, solve the l^{e−e₁} rest, and compose.
- #emph[Independent check:] the composite lands on E₂, and the degree-3⁷ isogeny is Kani-certified.

#strong[Baselines.] In `fixed_degree/baselines.rs`.

- #emph[Algorithms:] exhaustive search for prime d (Φ\_d roots, or kernels); MITM for d \= lᵉ; van Oorschot–Wiener claw search.
- #emph[Independent check:] planted paths are found and are cyclic. The counts are compared with theory in the MITM/vOW figure.

#strong[Quantum (Grover O(√d), Tani d^{1/3}).] Only a query-count model in `baselines::grover_tani_cost_model`, because no quantum computer is available.

=== √élu, the parts that were missing
#strong[§2 Strassen.] In `strassen.rs`.

- #emph[Algorithm:] ℓ! mod n for any n \< 2⁶², and deterministic factoring from gcd(n, ℓ! mod n).
- #emph[Independent check:] matches the naive product on 528 (ℓ, n) pairs; Wilson’s theorem at p ≈ 2²⁰; factoring matches `factor_u64`.

#strong[§3 Algorithm 1.] In `kernel/sqrt_mult.rs`.

- #emph[Algorithm:] the multiplicative h\_S over I + J ∪ K.
- #emph[Independent check:] equals the direct product over Zp, FpM\<2\> and FpM\<8\>; h\_S \= Φ\_n for 11 values of n.

#strong[§4.14 Irrational generators.] In `kernel/sqrt_velu_irr.rs`.

- #emph[Algorithm:] index systems made of Frobenius orbits; E\_J with coefficients in F\_q; norms taken down to F\_q.
- #emph[Independent check:] equals Montgomery Vélu over F\_{q^e} for e \= 2, 3, 4.

#strong[§6 CSURF.] In `path/csurf.rs`.

- #emph[Algorithm:] radical 2-steps (the formula is derived in the module) and odd steps by Vélu or √élu.
- #emph[Independent check:] every step matches Vélu and the modular polynomial Φ; the actions commute and invert at 62 bits, 111 bits and CSURF-512.

#strong[App. A.5.] In `kernel/sqrt_velu_mont.rs`.

- #emph[Algorithm:] fused Q(±1), Q(α), Q(1/α), and remainder trees with precomputed reciprocals.
- #emph[Independent check:] equals Vélu on CSIDH-512 and at ℓ \= 100003.

#strong[Tests.] `cargo test --release` in `isogeny_algos` passes 178 tests. The new test files are `tests/fixed_degree_pipeline.rs` (9 tests), `tests/fixed_degree_solvers.rs` (10), `tests/coppersmith.rs` (6) and `tests/sqrt_velu_gaps.rs` (15), plus 8 new unit tests.

== Measurements
The harness is `isogeny_algos/examples/fixed_degree_report.rs` (native Rust). #link("plot.py")[`plot.py`] only lays out the records.

=== End-to-end pipeline, p \= 1605631 \= 2¹⁵·49 − 1
These are 11 instances (figure "Wall time per step"). In every instance E₁ is the endpoint of an independent random 2-walk from 1728.

#emph[Planted prime d with √p \< d \< 2ᵃ \= 4096, Cornacchia (6 instances):]

- All 6 were found.
- The torsion representation was #strong[Kani-certified] in 6 of 6, and the pairing and Deuring-j checks passed in 6 of 6.
- d \= 1399 … 3671, so log\_p d \= 0.507 … 0.574.
- Total time 1.9 – 3.9 s per instance; 94–99 % of it is step 1 (83–87 % in the random-target instances below).

#emph[Random E₂, nothing planted, prime d ≈ 8–12p, so log\_p d ≈ 1.15–1.17 (3 instances):]

- All 3 were found: x \= \[111, 113, 14, 9\], \[−35, 249, 10, −4\], \[185, −65, 10, 8\].
- The norm equation took 226, 307 and 428 ms.
- The pairing and Deuring-j checks passed. Because d \> 2ᵃ there is no Kani certificate.

#emph[Coron’s bivariate method at this small p (2 instances):] #strong[no solution in the search box]. On the same instances Cornacchia found x \= \[0, −1, 0, 1\] and \[−1, 0, 1, −1\]. At this size the Gram entries are about 2√p ≈ 2¹¹, so W is too small for the lattice to certify any vector. See Limits.

The unit tests separately cover:

- planted smooth d (3·13·37 and 5·7·37 at p \= 1477631), producing a kernel representation;
- the hybrid method with d \= 3⁷;
- one Coron success at p \= 1605631 (`pipeline_coppersmith_bivariate`, where x \= \[−1, −1, 0, 0\]).

#figure([#image("figures/fig1_pipeline_steps.png")],
  caption: [
    Wall time per step per pipeline instance
  ]
)

#emph[Wall time of each step for each pipeline instance (run-20261016T0959Z, one run each). K \= Kani-certified torsion representation; p \= checked by pairing and Deuring j only.]

=== Steps 2–5 at cryptographic size
Each instance takes a random I₁ ⊂ O₀ and a planted J₀ of prime norm d. O₂ is the right order of a reduced ideal equivalent to I₁·J₀, so the Kirschmer–Voight ideal does not reveal d. There are two instances per setting (figure "Steps 2–5 wall time").

#strong[p \= 2⁶¹ − 1.] The table counts solved instances out of 2.

#align(center)[#table(
  columns: 5,
  align: (col, row) => (auto,auto,auto,auto,auto,).at(col),
  inset: 6pt,
  [d/√p], [log\_p d], [Cornacchia], [Coron], [trivariate],
  [≈ 2–9],
  [0.516–0.552],
  [—],
  [—],
  [2/6 overall: at d/√p \= 2 and 4; failed at 3, 5, 9, 9],
  [≈ 6],
  [0.542],
  [2/2],
  [2/2],
  [—],
  [≈ 20],
  [0.57],
  [2/2],
  [2/2],
  [—],
  [≈ 90],
  [0.606],
  [2/2],
  [2/2],
  [—],
  [≈ 380],
  [0.64],
  [2/2],
  [2/2],
  [—],
)
]

Wall time at log\_p d ≈ 0.64:

- Cornacchia 26.6 s and 71.2 s;
- Coron 6.5 s and 0.95 s.

Coron is faster than Cornacchia in 7 of the 8 (setting, repetition) pairs; each solver ran on its own random instances, so these are not paired on the same instance. The exception is rep 1 at d/√p ≈ 6 (0.019 s Cornacchia, 0.030 s Coron).

#strong[p \= 2¹²⁷ − 1, Coron only.] It solved 6 of 6 for d/√p \= 5–102 (log\_p d 0.518–0.553), taking 0.04–4.3 s.

Cornacchia is not run at 127 bits. Every guess needs an N of about 2²⁰⁰ factored, which Pollard rho cannot do; the slides’ L(1/3) cost assumes the number field sieve.

#figure([#image("figures/fig2_steps25_solvers.png")],
  caption: [
    Steps 2–5 wall time against log\_p d
  ]
)

#emph[Steps 2–5 wall time against log\_p d. Filled markers are 61-bit p, open diamonds 127-bit p. Dotted lines mark instances where the trivariate method found no solution.]

=== √élu §2: Strassen’s modular factorial
ℓ! mod (2⁶¹ − 1), mean over the repetitions stated in the records:

#align(center)[#table(
  columns: 4,
  align: (col, row) => (auto,auto,auto,auto,).at(col),
  inset: 6pt,
  [ℓ], [naive (ms)], [Strassen (ms)], [Strassen / naive],
  [10³],
  [0.044],
  [0.046],
  [1.05],
  [10⁴],
  [0.42],
  [0.43],
  [1.02],
  [10⁵],
  [4.06],
  [1.72],
  [0.42],
  [10⁶],
  [45.1],
  [11.7],
  [0.26],
  [4·10⁶],
  [164.8],
  [40.2],
  [0.24],
)
]

- Below b \= ⌊√ℓ⌋ \= 200, the code uses the naive product by design.
- The slope comes from Karatsuba multiplication: O(ℓ^0.79), not the paper’s Õ(√ℓ), which needs quasi-linear polynomial multiplication.

#figure([#image("figures/fig3_strassen.png")],
  caption: [
    Strassen vs naive modular factorial
  ]
)

=== Baselines for d \= 2ᵉ (p \= 4000039)
#strong[MITM.] Leaves equal 2·3·2^{e/2−1} exactly for e \= 12, 14, 16, 18 (figure "MITM and vOW counts").

#strong[vOW] (w \= 16) f-evaluations, against 2.5·√(|S|³/w):

#align(center)[#table(
  columns: 4,
  align: (col, row) => (auto,auto,auto,auto,).at(col),
  inset: 6pt,
  [e], [rep 0], [rep 1], [expected],
  [12],
  [1030],
  [1953],
  [1663],
  [14],
  [1056],
  [820],
  [4703],
  [16],
  [5933],
  [459],
  [13302],
  [18],
  [64744],
  [1016],
  [37624],
)
]

- The single-instance variance is large: the golden collision is found early or late.
- These spaces (|S| ≤ 1536) are too small to test the 1/√w scaling.

#figure([#image("figures/fig4_baselines.png")],
  caption: [
    MITM and vOW counts against theory
  ]
)

=== CSURF-512
One action with key entries in {−1, 0, 1} took 25.8 ms, and the inverse key returned to the starting curve.

== Limits and negative results
+ #strong[Prime size of the full pipeline.] The Deuring step (steps 1 and 6) is limited to p \< 2³¹ by `ext::Fp4`’s u128 arithmetic. It also needs a 256-smooth odd part of p² − 1 that is much larger than p; the primes were chosen by search. Steps 2–5 have no size limit.
+ #strong[Kani certificate.] It needs d \< 2ᵃ, where 2^(a+3) | p + 1. In this run that means d \< 4096 at p \= 1605631, i.e. log\_p d \< 0.58. For larger d the torsion images are checked only by the Weil pairing and the Deuring j.
+ #strong[Textbook auxiliaries do not work here (negative).] The usual auxiliary for the certificate is an endomorphism γ ∈ End(E₁) of degree 2ᵃ − d. In this regime it does not exist: 2ᵃ − d \< √p, and only integers have reduced norm that small in a maximal order. Robert’s dimension-4 integer matrix \[\[u, −w\], \[w, u\]\] was tried and failed: its first (2,2,2,2)-step is a product of two gluings, which `theta_g::chain_g` cannot start from. Integer α fails and Gaussian α works, checked on the crate’s own Kani-4 instances. The certificate actually used is a second quaternion-computed isogeny of degree 2ᵃ − d. Kani’s lemma needs only coprime degrees summing to 2ᵃ.
+ #strong[Coron needs large coefficients.] It failed in 2 of 2 pipeline instances at p ≈ 2²¹ where Cornacchia succeeded, while it worked in every 61- and 127-bit setting.
+ #strong[The trivariate method is weak and heuristic.] It succeeded in 2 of 6 settings at 61 bits (d/√p ≤ 4) and in the test at d/√p ≈ 4. The rigorous Bauer–Joux step is not implemented.
+ #strong[Cornacchia’s factoring.] It uses Pollard rho only, so N above about 2¹⁰⁰ with two large factors gets abandoned.
+ #strong[Strassen’s complexity.] It is Karatsuba-limited, as measured above.
+ #strong[Index systems.] §4.14 index systems for e ≥ 3, and Algorithm 1’s greedy index systems for unit groups with many factors, leave a large K. At the sizes tested, `SqrtVeluIrr` is slower than Vélu over the extension field, as the paper itself warns.
+ #strong[CSURF-512 prime.] It was reconstructed from the paper’s description. The paper’s brace says "6 consecutive primes" but the list has 5; readings that add 359 or 397 are not prime. The prime was not checked against the authors’ code.
+ #strong[Quantum algorithms.] Not simulated.
+ #strong[Lean formalisation.] There is none. The derived CSURF radical 2-step formula and the Kani-certificate argument are not formalised. The missing pieces are a Mathlib formalisation of Vélu’s formulas and of Kani’s lemma; the code’s correctness rests on the tests listed above. All other items are implementations of published algorithms.

== Reproduce
```bash
cd isogeny_algos
cargo test --release --test fixed_degree_pipeline --test fixed_degree_solvers --test coppersmith --test sqrt_velu_gaps
cargo run --release --example fixed_degree_report > ../research/fixed_degree_isogenies_20261016/data/records.jsonl
cd ../research/fixed_degree_isogenies_20261016 && python3 plot.py
FD_TRACE=1 cargo test --release --test fixed_degree_pipeline -- --nocapture   # per-step trace
```

== References
- B. Benčina, P. Kutas, S.-P. Merz, C. Petit, M. Stopar, C. Weitkämper, #emph[Finding isogenies of fixed degree between supersingular elliptic curves], Isogeny Club talk, 12 December 2023 (slides attached to the request).
- D. J. Bernstein, L. De Feo, A. Leroux, B. Smith, #emph[Faster computation of isogenies of large prime degree], 2020 (PDF attached to the request).
- J.-S. Coron, #emph[Finding small roots of bivariate integer polynomial equations: a direct approach], CRYPTO 2007.
- A. Bauer, A. Joux, #emph[Toward a rigorous variation of Coppersmith’s algorithm on three variables], EUROCRYPT 2007.
- M. Kirschmer, J. Voight, #emph[Algorithmic enumeration of ideal classes for quaternion orders], SIAM J. Comput. 2010.
- K. Eisenträger, S. Hallgren, K. Lauter, T. Morrison, C. Petit, #emph[Supersingular isogeny graphs and endomorphism rings: reductions and solutions], EUROCRYPT 2018.
- D. Kohel, K. Lauter, C. Petit, J.-P. Tignol, #emph[On the quaternion ℓ-isogeny path problem], LMS J. Comput. Math. 2014.
- D. Robert, #emph[Evaluating isogenies in polylogarithmic time], 2022.
- W. Castryck, T. Decru, #emph[CSIDH on the surface], PQCrypto 2020.
- G. Cornacchia, #emph[Su di un metodo per la risoluzione in numeri interi dell’equazione ∑ c\_h x^{n−h} y^h \= P], 1908.

These references are cited from memory, except the two attached sources, which were read in this session.
