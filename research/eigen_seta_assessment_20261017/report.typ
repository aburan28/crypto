#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#set heading(numbering: "1.")
#align(center)[#text(size: 15pt, weight: "bold")[Eigen-Séta v2 assessment, 2026-10-17]\
Authors' own analysis, not an external review. Code: `eigenx` commit `94fe2d9`. Paper: `crypto:papers/eigen-seta/` revision 2.]

= Verdict
No attack below the stated costs was found. The binding attack is still mask guessing: $approx 2^lambda$ classical chains, and $2^(lambda\/2)$ Grover queries on L1 and L3, $>= 2^lambda$ queries on L5, Q3 and Q5. Risks: (1) an effective Z[f ω]-orientation of $E_A$, which would give subexponential quantum key recovery at *every* set; (2) a mask-free dimension-4/8 embedding; (3) the decisional assumption.

= Findings
#table(columns: (auto, 1fr, auto), stroke: 0.4pt, inset: 4pt,
  [*\#*], [*Finding*], [*Kind*],
  [F1], [Candidate $f$ testable from the public key alone: Kani chain on $E_A times E_A$ of $e approx 2 log_2 f + 3$ steps vs $E approx 4 log_2 f$ for decapsulation. Degree-search gate columns −1.9 bits; no tier changes.], [proved, measured],
  [F2], [Revision 1's NIST categories 1/3/5 for L1/L3/L5 hold only for gate count under MAXDEPTH; L1/L3 have $2^64.5$ / $2^96.5$ quantum queries. Claim withdrawn in unqualified form.], [derived],
  [F3], [Hidden orientation is the largest structural risk (Kuperberg if effective); Q sets enlarge $d$, not $f$.], [literature],
  [F4], [`Fp::mul2/mul3` wrong for non-friendly moduli; `tests/field.rs::f1` red on main; no parameter set affected. Fixed.], [measured],
  [F5], [SECURITY.md placed the ω-eigenbasis on $E_D[3^F]$ (impossible: 3 inert in $QQ(sqrt(-7))$). Corrected.], [derived],
  [F6], [Model filter fraction $1\/ln 2^E$ mis-justified; measured ρ is 1.6–2.3 bits higher (model conservative).], [derived],
  [F7], [eigenx CI red on main at fmt/clippy. Fixed.], [measured],
)

= Lemma 1 (public-key test)
With $psi_A = 2 tau - t f$ (trace 0, degree $7 f^2$), odd $u$, and $2^e - 7 u^2 f'^2 = c^2 + c'^2$: the kernel $⟨(c' P', -lambda_S P'), (c' Q', -lambda_T Q')⟩$ is isotropic for every candidate $f'$, since $c'^2 + lambda_S lambda_T equiv c^2 + c'^2 + 7u^2 f'^2 = 2^e$. For $f' = f$ it is Kani's kernel for $gamma = u psi_A + c$, which splits into $E_A times E_A$. Lean (core, v4.34.1): `pkTest_kernel_isotropic`, `pkTest_kani_degree` pass the proof gate; the bridge from norm to isogeny degree is not formalised.

= Measurement: run-20261017T1623Z-pk-degree-test
Shared 4-vCPU Xeon @ 2.10 GHz VM, rustc 1.97.0, no pinning. Medians [ranges].
#table(columns: 8, stroke: 0.4pt, inset: 3pt,
  [set], [$E$], [$e$ true], [$e_min$], [excess], [wrong (acc.)], [test s], [decaps s],
  [L1], [1076], [569], [560], [9.9], [105 (0)], [0.79 [0.53–0.91] n=7], [1.28 [1.11–1.57]],
  [L3], [1594], [832], [823], [14.7], [35 (0)], [2.72 [2.64–2.78] n=3], [5.54 [5.49–5.75]],
  [L5], [2106], [1087], [1076], [14.0], [18 (0)], [6.84 [6.66–7.03] n=3], [13.34 [11.54–14.26]],
)
#figure(image("figures/pk_test.svg", width: 95%), caption: [Chain length and wall time of the public-key test vs decapsulation; measured, run run-20261017T1623Z-pk-degree-test.])
#figure(image("figures/attack_margins.svg", width: 85%), caption: [Margins over λ by tier, derived by the eigenx cost model (records/security_estimate.txt).])
#figure(image("figures/attack_map.svg", width: 100%), caption: [Structure of the scheme and the attack routes considered (diagram, not a measurement).])

= Unresolved
Risks 1–3 above (eigenx LITERATURE_DRILL L2, L6, L7, L8); external review of the $2^8$ Grover-oracle constant; constant-time evidence on the theta build. Records: `records/`. Markdown twin: `REPORT.md`.
