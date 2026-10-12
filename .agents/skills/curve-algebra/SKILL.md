---
name: curve-algebra
description: Produce source-linked technical explanations, algebraic diagrams, and illustrated PDFs for elliptic-curve structure. Use when documenting CM discriminants, endomorphisms, GLV/GLS decomposition, Frobenius, automorphism actions, or the connection between formulas and algorithmic costs. Preserve verification and evidence gates; documentation does not establish a new speedup.
---

# Curve algebra

Build a report that connects each algebraic object to its algorithmic role.

## Workflow

1. Read the visible question and applicable repository instructions. Identify the missing bridge: curve data to CM order, order to map, map to subgroup eigenvalue, eigenvalue to lattice, or action to orbit quotient. Proceed with the established examples when context is sufficient.
2. Read [math-contract.md](references/math-contract.md) for notation, scope, exact fixtures and primary sources. Verify unfamiliar or source-specific claims against primary publications. Distinguish the equation discriminant, Frobenius discriminant, order discriminant and map-polynomial discriminant.
3. Read [report-layout.md](references/report-layout.md). Adapt its page plan and editable diagrams to the question. Explain one concept per diagram; place equations beside the step they justify. Show a positive example and a failure case. Label synthetic examples and distinguish them from standardized curves or prior benchmarks.
4. Run the native fixture checker before citing its examples. Compile into scratch, never into this skill directory:

   ```sh
   c++ -std=c++17 -O2 -Wall -Wextra -Werror scripts/verify_examples.cpp -o /tmp/curve-algebra-check
   /tmp/curve-algebra-check
   ```

   Interpret the output as an arithmetic verification receipt, not a performance benchmark. For different examples, use the project's native arithmetic and verifier; do not introduce Python curve arithmetic or research harnesses.
5. Produce the requested diagrams and PDF. Prefer vector shapes and typeset equations using an available document renderer. Use the explicitly selected PDF plugin if it exposes relevant callable capabilities; otherwise explain its absence briefly and use an available local renderer. Restrict any Python used to document layout, formula rendering and PDF inspection.
6. Add primary-source links and short questions with answers. Keep geometric degree, map evaluation cost, scalar-multiplication improvement, ideal rho iteration reduction and measured ECDLP cost as different quantities. Mark unmeasured gains as unknown.
7. Render every page for visual inspection. Check legibility, symbols, arrows, clipping, overlaps, citations and page breaks. Confirm the PDF opens. Reuse completed arithmetic validation; do not claim timings from fixture checks.
8. Save the user-facing deliverable using the active file-saving workflow; replace an existing file's identity when updating it. For a repository task, preserve editable source and follow its PR, evidence and identity rules. Return the PDF link with a suggested starting page. Do not return a packaged skill unless requested.

## Scope and claim gates

- State ordinary-curve and characteristic assumptions where the equations require them.
- Establish the actual order before treating its elements as maps on the named curve.
- Verify field of definition, subgroup preservation and φ(P) = [λ]P. A polynomial root alone is a candidate eigenvalue.
- Explain that a balanced short kernel-lattice basis supports GLV bounds; its determinant alone does not.
- Describe rho folding through an effective cheap action, mostly full-sized orbits and a compatible walk. Present √m as an ideal birthday iteration factor, with canonicalization and cycle costs still to price.
- For higher dimensions, name the compatible maps and the applicability of the published construction. Four coefficients alone do not prove an r^(1/4) bound or a fourfold wall-time gain.
- For index calculus, retain relative phases and compatibility costs. Orbit compression alone does not establish cheaper relations.
- When a technical task turns into a substantive new search, use the repository's research-visuals workflow and full-cost evidence gates. Keep algebraic fixtures separate from benchmark results and canonical scoreboards.
