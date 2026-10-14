PRAGMA foreign_keys = ON;
CREATE TABLE source (
  source_id TEXT PRIMARY KEY,
  title TEXT NOT NULL,
  url TEXT NOT NULL,
  source_scope TEXT NOT NULL
) STRICT;
INSERT INTO source VALUES (
  'IMC2015',
  'A classification of elliptic curves with respect to the GHS attack in odd characteristic',
  'https://eprint.iacr.org/2015/805.pdf',
  'Section 5 numbered cases; odd characteristic, elliptic source, prime extension degree 2, 3, 5, or 7, (2,...,2)-covering'
);
CREATE TABLE published_case (
  case_no INTEGER PRIMARY KEY CHECK (case_no BETWEEN 1 AND 36),
  d INTEGER NOT NULL CHECK (d IN (2, 3, 5, 7)),
  n INTEGER NOT NULL CHECK (n BETWEEN 2 AND d),
  e INTEGER NOT NULL CHECK (e >= 0),
  cover_genus INTEGER NOT NULL CHECK (cover_genus = d + e),
  coefficient TEXT NOT NULL CHECK (coefficient IN ('one', 'eta')),
  printed_page INTEGER NOT NULL CHECK (printed_page IN (17, 18)),
  source_id TEXT NOT NULL REFERENCES source(source_id)
) STRICT;
CREATE VIEW case_profile AS SELECT
  case_no, d, n, e, cover_genus, coefficient, printed_page, source_id,
  CASE WHEN e = 0 THEN 'genus_minimal' ELSE 'higher_genus' END AS covering_profile,
  CASE WHEN coefficient = 'eta' THEN 'square_or_nonsquare' ELSE 'square_only' END AS permitted_square_class
FROM published_case;
CREATE TABLE minimum_case_pattern (
  case_no INTEGER PRIMARY KEY REFERENCES published_case(case_no),
  hd_pattern TEXT NOT NULL,
  h1_degrees TEXT NOT NULL,
  root_scope TEXT NOT NULL
) STRICT;
CREATE VIEW minimum_case_recipe AS SELECT
  p.case_no, p.d, p.n, p.cover_genus, p.permitted_square_class,
  m.hd_pattern, m.h1_degrees, m.root_scope
FROM case_profile AS p JOIN minimum_case_pattern AS m USING (case_no)
WHERE p.e = 0;
