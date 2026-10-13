# Hu, ePrint 2025/015: access record and implementability

**Paper.** Jianjun Hu, *A New Method for Solving Discrete Logarithm Based on
Index Calculus*, IACR ePrint 2025/015 (received 2025-01-04, revised
2025-09-02, preprint).
**Companion.** Jianjun Hu, *Implementation analysis of index calculus method on
elliptic curves over prime finite fields*, IACR ePrint 2024/1923 (2024-11-27;
also EasyChair Preprint 15613, 2024-12-20).

## What was accessible (2026-10-11)

| Source | Result |
|---|---|
| `curl https://eprint.iacr.org/2025/015.pdf` | HTML page titled "Just a moment..." (Cloudflare bot challenge), not a PDF. Not bypassed. |
| WebFetch `https://eprint.iacr.org/2025/015.pdf` | HTTP 403. Not bypassed. |
| WebFetch `https://eprint.iacr.org/2025/015` (abstract page) | Served: title, author, dates, abstract. |
| `curl https://eprint.iacr.org/2024/1923.pdf` | Same Cloudflare challenge page. |
| WebFetch `https://eprint.iacr.org/2024/1923` (abstract page) | Served. |
| Web search for other copies of 2025/015 (EasyChair, arXiv, mirrors) | None found; the EasyChair preprint found is the 2024/1923 paper. |
| Wayback Machine availability API | HTTP 429 (rate limited); no copy obtained. |

The full text of 2025/015 was therefore **not obtained**. This is an access
limit, not a judgement about the paper. The brief said not to bypass the
challenge, and it was not bypassed.

## What the abstract says

- The subject is the discrete logarithm problem **in prime finite fields**, by
  index calculus (IC). It is not an elliptic-curve method. The literature
  search filed it under "Claimed ECDLP attacks"; from the abstract, that filing
  is wrong. The paper is about `F_p^*`.
- It says an earlier method, "IICA", studied only the case "integer up" and
  not "integer down". It says the down direction needs modular multiplication.
  The abstract never expands "IICA" and never defines either direction.
- The claims are that the success probability of solving the discrete
  logarithm nearly doubles, that fewer transformations are needed, and that
  time complexity drops. No sizes, no parameters and no asymptotic exponent
  are given. The abstract says this poses a "serious challenge" to
  prime-field cryptosystems.

## Why nothing could be implemented or measured

1. Without the full text, the algorithm is unknown. "IICA", "integer up",
   "integer down" and the "transformation" being counted are undefined in
   every source that could be reached. Any implementation would be a guess at
   what the paper means. It would not reproduce the paper.
2. Even taken at face value, the claim is a constant factor (about 2× in
   success probability) inside `F_p^*` index calculus. That is outside this
   repository's ECDLP scope. It says nothing about elliptic curves over prime
   fields.
3. The same author's ePrint 2024/1923 (abstract) states that, as far as
   current research knows, no subexponential method for the ECDLP over prime
   fields exists. Its "pseudo-index calculus" builds on Petit et al. (2016) and
   Momonari–Kudo et al. (2018) and claims correctness only, not an advance in
   complexity.

**Status: blocked (full text inaccessible without bypassing the bot
challenge).** The original requirement, "implement and measure if it specifies
an algorithm", stays unresolved. The finding that it targets `F_p^*` and not
the ECDLP comes from the abstract alone.
