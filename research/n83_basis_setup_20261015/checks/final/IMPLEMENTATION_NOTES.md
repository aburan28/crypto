# Retained implementation and publication failures

No scientific retry was performed. First new C++ aggregation build required a
std::string-to-c_str conversion for the inherited require(bool,const char*) API;
the test used Parser construction/parse instead of inherited aggregate/value.
Both compile issues were corrected only in new postprocessing source. Successful
coverage, raw-SHA and recursive A/A tamper tests are recorded. No frozen v1
source or binary changed.

Initial site/browser validation lacked sparse inputs. Source paths were
materialized, including the browser's full committed discovery surface, before
regeneration. A sparse browser rendering was rejected because it would lose old
sessions/candidates. Final counts equal historical counts. New scoreboard panel
needed its panel-summary/index entry and plain<tr><td> plotted-row convention;
presentation was corrected without changing tests. All48 final site tests pass.
Earlier failure logs remain.

Broad curve-name linter attempts encountered unrelated off-disk documents; the
unchanged linter was run on every changed/new Markdown/HTML document instead.
That targeted check passes; no linter guard/allowlist was changed. Normal full-
checkout CI retains its global check. Stage tests/library/source seals and PDF
preservation results are separate, not inferred from that targeted check.
