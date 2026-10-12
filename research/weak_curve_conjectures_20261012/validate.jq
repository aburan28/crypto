def nonempty: type == "string" and length > 0;
if .schema != "weak_curve.conjectures/v1" then error("schema") else . end
| . as $root
| if (.entries | length) != 8 then error("entry count") else . end
| if ([.entries[].id] | unique | length) != 8 then error("duplicate id") else . end
| if ([.entries[] | select(.status == "refuted")] | length) != 1 then error("refuted count") else . end
| if ([.entries[] | select(.status == "open")] | length) != 7 then error("open count") else . end
| if all(.entries[];
    (.id | startswith("WCQ1:")) and
    (.title | nonempty) and (.statement | nonempty) and
    (.domain | nonempty) and (.novelty | nonempty) and
    (.falsifier | nonempty) and (.next_experiment | nonempty) and
    (.lean_status | nonempty) and
    all(.references[]; . as $id | any($root.references[]; .id == $id)))
  then "CONJECTURE_LEDGER_VALID|entries=8|open=7|refuted=1"
  else error("incomplete hypothesis or unknown source") end
