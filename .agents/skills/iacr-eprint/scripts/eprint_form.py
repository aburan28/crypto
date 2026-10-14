#!/usr/bin/env python3
"""Prepare an IACR ePrint submission from a LaTeX paper directory.

Writes two files next to the paper (or into --out):
  eprint-abstract.txt  the abstract as the form accepts it
  EPRINT.md            every form field with its value, plus the checks

and exits non-zero if a hard ePrint rule fails (anonymous first page,
no email in the PDF, bad abstract, bad keywords, wrong paper size).

  python3 eprint_form.py paper/ecc2k130 \
      --category "Attacks and cryptanalysis" \
      --keywords "ECC2K-130, Pollard rho, Koblitz curves" \
      --license "CC BY"

Rules encoded here come from https://eprint.iacr.org/submit and
https://eprint.iacr.org/operations.html (read 2026-10-13); re-read them
if the form changes.
"""
import argparse, pathlib, re, subprocess, sys, unicodedata

CATEGORIES = ["Applications", "Cryptographic protocols", "Foundations",
              "Implementation", "Secret-key cryptography",
              "Public-key cryptography", "Attacks and cryptanalysis"]
LICENSES = ["CC BY", "CC BY-SA", "CC BY-NC", "CC BY-NC-SA", "CC BY-NC-ND", "CC0"]
# Text-mode commands the form rejects (math-mode LaTeX inside $...$ is fine).
# Formatting ones are unwrapped first (UNWRAP); anything left here is a FAIL.
TEXT_CMDS = r"\\(emph|textbf|textit|texttt|textsc|textrm|textsf|mbox|cite|citep|citet|ref|eqref|label|footnote|url|href|art|section|item)\b"
# Paper macros worth expanding so MathJax can render them.
MACROS = {r"\FF": r"\mathbb{F}", r"\ZZ": r"\mathbb{Z}", r"\QQ": r"\mathbb{Q}",
          r"\NN": r"\mathbb{N}", r"\RR": r"\mathbb{R}", r"\CC": r"\mathbb{C}"}

def read_tex(d):
    for name in sorted(d.glob("*.tex")):
        s = name.read_text()
        if r"\documentclass" in s:
            return name, s
    sys.exit(f"no main .tex (with \\documentclass) in {d}")

def braced(s, cmd):
    m = re.search(r"\\" + cmd + r"\s*\{", s)
    if not m:
        return None
    i, depth = m.end(), 1
    while depth:
        depth += {"{": 1, "}": -1}.get(s[i], 0)
        i += 1
    return s[m.end():i - 1]

def strip_comments(s):
    return re.sub(r"(?<!\\)%.*", "", s)

def abstract_source(d, tex):
    m = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    if not m:
        sys.exit("no abstract environment")
    body = m.group(1)
    inp = re.search(r"\\input\{([^}]+)\}", body)
    if inp:
        body = (d / inp.group(1)).read_text()
    return strip_comments(body)

# Formatting-only text commands: the form rejects them, but their content
# reads fine as plain text, so unwrap them (outside math only).
UNWRAP = r"\\(emph|textbf|textit|texttt|textsc|textrm|textsf|mbox)\{([^{}]*)\}"

def unwrap_outside_math(s):
    parts = re.split(r"(\$[^$]*\$)", s)
    return "".join(p if p.startswith("$") else re.sub(UNWRAP, r"\2", p) for p in parts)

def to_form_text(s):
    s = unwrap_outside_math(s)
    for k, v in MACROS.items():
        s = re.sub(re.escape(k) + r"(?![A-Za-z])", lambda _: v, s)
    s = s.replace("~", " ").replace("{,}", ",").replace("---", "—").replace("--", "–")
    s = s.replace("``", "“").replace("''", "”").replace(r"\%", "%")
    s = re.sub(r"\s+", " ", s).strip()
    return s

def plain_title(t):
    t = strip_comments(t).replace("\\\\", " ").replace("~", " ").replace("--", "–")
    return re.sub(r"\s+", " ", t).strip()

def authors_plain(a):
    a = strip_comments(a)
    a = re.sub(r"\\texttt\{([^}]*)\}", r"\1", a)
    parts = re.split(r"\\\\", a)
    names = [n.strip() for n in re.split(r"\\and", parts[0]) if n.strip()]
    rest = [re.sub(r"\s+", " ", p).strip() for p in parts[1:] if p.strip()]
    return names, rest

def pdf_text(pdf, page=1):
    return subprocess.run(["pdftotext", "-f", str(page), "-l", str(page), str(pdf), "-"],
                          capture_output=True, text=True).stdout

def pdf_info(pdf):
    out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True).stdout
    return dict(l.split(":", 1) for l in out.splitlines() if ":" in l)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paper_dir")
    ap.add_argument("--category", required=True, choices=CATEGORIES)
    ap.add_argument("--keywords", required=True, help="comma-separated")
    ap.add_argument("--license", default="CC BY", choices=LICENSES)
    ap.add_argument("--pdf", help="defaults to <main>.pdf")
    ap.add_argument("--publication", default="Published nowhere else")
    ap.add_argument("--out", help="output dir (default: paper dir)")
    a = ap.parse_args()

    d = pathlib.Path(a.paper_dir)
    out = pathlib.Path(a.out or d)
    main_tex, tex = read_tex(d)
    pdf = pathlib.Path(a.pdf) if a.pdf else main_tex.with_suffix(".pdf")
    errors, warnings = [], []

    title = plain_title(braced(tex, "title") or "")
    if not title:
        errors.append("no \\title")
    if re.search(r"<[a-zA-Z/]|&", title):
        warnings.append("title contains '&' or HTML; the form discourages both")

    names, affil = authors_plain(braced(tex, "author") or "")
    if not names or any(n.lower() == "anonymous" for n in names):
        errors.append("\\author is Anonymous: ePrint takes no anonymous papers (restore the author block)")
    if any("\\" in n or "{" in n for n in names):
        errors.append("author names must be plain UTF-8, not TeX: " + ", ".join(names))

    abstract = to_form_text(abstract_source(d, tex))
    bad = sorted(set(m.group(0) for m in re.finditer(TEXT_CMDS, abstract)))
    if bad:
        errors.append("abstract uses text-mode LaTeX the form rejects: " + " ".join(bad))
    if re.search(r"<[a-zA-Z/]", abstract):
        errors.append("abstract contains HTML tags")
    if abstract.count("$") % 2:
        errors.append("abstract has an unbalanced $")
    if len(abstract) < 64:
        errors.append(f"abstract is {len(abstract)} chars; the form needs at least 64")
    nonprint = [c for c in abstract if unicodedata.category(c) in ("Cc", "Cf", "Co", "Cn")]
    if nonprint:
        errors.append(f"abstract has {len(nonprint)} non-printable characters")
    leftover = sorted(set(re.findall(r"\\[A-Za-z]+", re.sub(r"\$[^$]*\$", "", abstract))))
    if leftover:
        warnings.append("commands outside math (check they render): " + " ".join(leftover))

    kws = [k.strip() for k in a.keywords.split(",") if k.strip()]
    for k in kws:
        if len(k) > 40:
            errors.append(f"keyword over 40 chars: {k!r}")
        if "\\" in k or "$" in k:
            errors.append(f"keyword contains LaTeX: {k!r}")

    if not pdf.exists():
        errors.append(f"{pdf} missing: build it first")
    else:
        p1 = pdf_text(pdf)
        info = pdf_info(pdf)
        if re.search(r"\banonymous\b", p1, re.I):
            errors.append("page 1 of the PDF says Anonymous")
        if not re.search(r"[\w.+-]+@[\w-]+\.[\w.-]+", p1):
            errors.append("no email address on page 1 of the PDF (the form requires author emails in the PDF)")
        for n in names:
            if n and n not in re.sub(r"\s+", " ", p1):
                warnings.append(f"author {n!r} not found verbatim on page 1")
        if info.get("Author", "").strip().lower() in ("", "anonymous"):
            warnings.append("PDF metadata Author is empty or Anonymous")
        m = re.search(r"([\d.]+) x ([\d.]+) pts", info.get("Page size", ""))
        if m:
            w, h = float(m.group(1)), float(m.group(2))
            if not any(abs(w - W) < 3 and abs(h - H) < 3 for W, H in [(595.28, 841.89), (612, 792)]):
                errors.append(f"page size {w}x{h} pt is neither A4 nor US Letter")

    emails = re.findall(r"[\w.+-]+@[\w-]+\.[\w.-]+", " ".join(affil))
    if emails and len(emails) == len(names):
        author_cell = "; ".join(f"{n} <{e}>" for n, e in zip(names, emails))
    else:
        author_cell = "; ".join(names)
        if not emails:
            warnings.append("no email in the \\author block; the form needs at least one")
    (out / "eprint-abstract.txt").write_text(abstract + "\n")
    md = f"""# IACR ePrint submission: {title}

Generated by `.agents/skills/iacr-eprint/scripts/eprint_form.py`; re-run it
after any change to the title, authors, abstract or PDF. Submit at
https://eprint.iacr.org/submit. **Only an author submits**: the license is
irrevocable, an accepted paper can never be fully removed, and every named
author must agree.

| field | value |
|---|---|
| Title | {title} |
| Authors (in order, UTF-8) | {author_cell} |
| Affiliation / contact (on page 1) | {"; ".join(affil)} |
| Abstract | paste `eprint-abstract.txt` ({len(abstract)} chars; `$…$` math renders with MathJax) |
| Category | {a.category} |
| Keywords | {", ".join(kws)} |
| Publication information | {a.publication} |
| Content | `{pdf.name}` ({pdf_info(pdf).get("Pages", "?").strip() if pdf.exists() else "?"} pages) |
| Additional notes (public) | empty for a first version; on a revision, what changed |
| Message to editors (private) | empty unless something needs explaining |
| License | {a.license} |
| Confirmation | tick after reading the acceptance, license and withdrawal conditions |

Enter each author with an email where known (at least one is required);
ORCID and ROR identifiers are optional.

## Checks

{"".join(f"- FAIL: {e}" + chr(10) for e in errors) or "- all hard rules pass" + chr(10)}{"".join(f"- warn: {w}" + chr(10) for w in warnings)}
## After acceptance

Keep the password from the acceptance email. Revisions go through
https://eprint.iacr.org/revise on the same number, never as a new entry.
Record the ePrint number in the paper's README.
"""
    (out / "EPRINT.md").write_text(md)
    for e in errors:
        print("FAIL:", e)
    for w in warnings:
        print("warn:", w)
    print(f"wrote {out/'eprint-abstract.txt'} and {out/'EPRINT.md'}")
    sys.exit(1 if errors else 0)

if __name__ == "__main__":
    main()
