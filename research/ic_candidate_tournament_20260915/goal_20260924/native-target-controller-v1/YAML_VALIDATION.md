# Workflow parser compatibility observation

The first local YAML check used Ruby's `YAML.load_file(ARGV[0], aliases: true)`.
It exited 1 before parsing because this host's system Ruby 2.6 does not accept
that keyword. Original error:

```text
/System/Library/Frameworks/Ruby.framework/Versions/2.6/usr/lib/ruby/2.6.0/psych.rb:576:in `load_file': unknown keyword: aliases (ArgumentError)
        from -e:1:in `<main>'
```

The compatible check uses `YAML.load_file(ARGV[0])`; its output is retained in
`yaml-v2.log`. The workflow source was not changed to address this local parser
API mismatch. Both invocations used the native busy wrapper. This is workflow
syntax evidence only, not a research execution or GitHub CI pass.
