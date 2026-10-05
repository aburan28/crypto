# Original target freeze failure retained

The first target freeze was invoked from source commit `98935a28f` with the
frozen `config.json`, `preparation-binding.json` and `host-context.json` in this
directory. It exited 1 before creating the target capsule or invoking a target
worker. Its original error was:

```text
icprog: invalid type: floating point `0.007025581`, expected JSON with unique keys and integer numbers at line 24 column 37
```

The field was `log_report.collection_seconds` in the accepted F5 ordinary
preparation producer. It is a diagnostic float outside `mathematical_input`;
the target checker had incorrectly applied its integer-only target-result
parser to the whole ordinary producer. The subsequent source correction uses a
separate duplicate-key-safe parser that permits finite diagnostic floats for
that one source-bound preparation record. Target result records remain
integer-only. This failed build is not a scientific execution or a registration
consumption, and no target evidence or online time can be derived from it.
