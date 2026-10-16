# delphi has moved to Conductor

`delphi` delegates prompts to frontier models through provider batch APIs (Anthropic
Message Batches, OpenAI Batch): half price, answered minutes to hours later. It now lives in
[aburan28/conductor](https://github.com/aburan28/conductor) (added in
[aburan28/conductor#60](https://github.com/aburan28/conductor/pull/60); until it merges, use
that branch, `claude/delphi-delegation`):

- **Hosted by Conductor** — `conductor delphi ask --task T-42 "..."`, keys from the project
  vault, spend in `conductor usage`, a `delphi.finished` event on the linked task. See
  `docs/DELPHI.md` there.
- **Standalone binary** — no control plane, keys from the environment:

  ```sh
  go install github.com/aburan28/conductor/cmd/delphi@main
  export ANTHROPIC_API_KEY=...            # and/or OPENAI_API_KEY + DELPHI_OPENAI_MODEL
  delphi run -prompt "Say hello" -out runs/smoke
  delphi serve &                           # the loop; then: delphi ask / delphi get -wait
  ```

Default model: `claude-mythos-5-1`. Not yet verified against the live provider APIs.
