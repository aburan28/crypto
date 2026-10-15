package main

// Request files: parsing, defaults and local validation. Everything here runs
// without network or credentials. The checks encode the constraints of
// claude-mythos-5-1 so a bad row fails before a billable batch is created:
// invalid rows otherwise surface only after the asynchronous round trip.

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"regexp"
	"sort"
	"strconv"
)

const (
	defaultModel = "claude-mythos-5-1"
	// Thinking is always on and counts toward max_tokens, so leave room for it.
	// Batches are asynchronous: no HTTP timeout pushes this down.
	defaultMaxTokens = 32000
	maxOutputTokens  = 128000
)

// Per-batch limits; variables so tests can shrink them. The API allows 100,000
// requests and 256 MB per batch; stay under the byte limit with headroom.
var (
	maxRequestsPerBatch = 100000
	maxBatchBytes       = 200 << 20
)

var (
	efforts      = []string{"low", "medium", "high", "xhigh", "max"}
	customIDRe   = regexp.MustCompile(`^[A-Za-z0-9_-]{1,64}$`)
	rowKeys      = map[string]bool{"custom_id": true, "prompt": true, "messages": true, "system": true, "max_tokens": true, "effort": true, "params": true}
	forbiddenKey = map[string]string{
		"temperature": "sampling parameters are not accepted by Claude Mythos 5.1",
		"top_p":       "sampling parameters are not accepted by Claude Mythos 5.1",
		"top_k":       "sampling parameters are not accepted by Claude Mythos 5.1",
		"stream":      "batch requests are non-streaming",
		"model":       "set the model with -model / MYTHOS_MODEL, not per row",
		"max_tokens":  "use the row-level max_tokens key",
		"messages":    "use the row-level prompt or messages key",
		"system":      "use the row-level system key",
	}
)

// Defaults apply to rows that do not set the field themselves.
type Defaults struct {
	Model     string
	MaxTokens int
	Effort    string
	System    string
}

func defaultsFromEnv() (Defaults, error) {
	d := Defaults{Model: defaultModel, MaxTokens: defaultMaxTokens}
	if v := os.Getenv("MYTHOS_MODEL"); v != "" {
		d.Model = v
	}
	if v := os.Getenv("MYTHOS_MAX_TOKENS"); v != "" {
		n, err := strconv.Atoi(v)
		if err != nil {
			return d, fmt.Errorf("MYTHOS_MAX_TOKENS=%q is not an integer", v)
		}
		d.MaxTokens = n
	}
	return d, nil
}

func (d Defaults) check() []string {
	var p []string
	if d.MaxTokens < 1 || d.MaxTokens > maxOutputTokens {
		p = append(p, fmt.Sprintf("max_tokens must be in 1..%d", maxOutputTokens))
	}
	if d.Effort != "" && !contains(efforts, d.Effort) {
		p = append(p, fmt.Sprintf("effort must be one of %v, got %q", efforts, d.Effort))
	}
	return p
}

// Request is one entry of requests.jsonl: the API request, ready to submit.
type Request struct {
	CustomID string         `json:"custom_id"`
	Params   map[string]any `json:"params"`
}

type requestError struct{ Problems []string }

func (e *requestError) Error() string {
	return fmt.Sprintf("%d problem(s): %s", len(e.Problems), e.Problems[0])
}

func contains(xs []string, s string) bool {
	for _, x := range xs {
		if x == s {
			return true
		}
	}
	return false
}

func textBlocks(s string) []any {
	return []any{map[string]any{"type": "text", "text": s}}
}

func checkParams(extra map[string]any, where string, add func(string, ...any)) {
	keys := make([]string, 0, len(extra))
	for k := range extra {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		if why, bad := forbiddenKey[k]; bad {
			add("%s: params.%s: %s", where, k, why)
		}
	}
	if t, ok := extra["thinking"]; ok {
		m, _ := t.(map[string]any)
		if m == nil || m["type"] != "adaptive" {
			add("%s: params.thinking must be omitted or {\"type\":\"adaptive\"} (always-on thinking; disabled/enabled+budget_tokens are a 400)", where)
		}
	}
	if tc, ok := extra["tool_choice"].(map[string]any); ok {
		if t := tc["type"]; t == "any" || t == "tool" {
			add("%s: params.tool_choice type %q is a 400 on Claude Mythos 5.1; use \"auto\" (or \"none\") and name the tool in the prompt", where, t)
		}
	}
}

// normalizeRow validates one input row and returns the API request.
func normalizeRow(row map[string]any, d Defaults, where string) (Request, []string) {
	var problems []string
	add := func(f string, a ...any) { problems = append(problems, fmt.Sprintf(f, a...)) }

	var unknown []string
	for k := range row {
		if !rowKeys[k] {
			unknown = append(unknown, k)
		}
	}
	if len(unknown) > 0 {
		sort.Strings(unknown)
		add("%s: unknown keys %v (allowed: custom_id prompt messages system max_tokens effort params)", where, unknown)
	}

	cid, _ := row["custom_id"].(string)
	if !customIDRe.MatchString(cid) {
		add("%s: custom_id must match [A-Za-z0-9_-]{1,64}, got %v", where, row["custom_id"])
	}

	_, hasPrompt := row["prompt"]
	_, hasMessages := row["messages"]
	var messages []any
	switch {
	case hasPrompt == hasMessages:
		add("%s: give exactly one of 'prompt' or 'messages'", where)
	case hasPrompt:
		if s, ok := row["prompt"].(string); !ok || len(bytes.TrimSpace([]byte(s))) == 0 {
			add("%s: prompt must be a non-empty string", where)
		} else {
			messages = []any{map[string]any{"role": "user", "content": textBlocks(s)}}
		}
	default:
		raw, _ := row["messages"].([]any)
		if len(raw) == 0 {
			add("%s: messages must be a non-empty list", where)
			break
		}
		for i, m := range raw {
			mm, _ := m.(map[string]any)
			role, _ := mm["role"].(string)
			if mm == nil || (role != "user" && role != "assistant" && role != "system") {
				add("%s: messages[%d] needs role user, assistant or system", where, i)
				continue
			}
			out := map[string]any{"role": role}
			switch c := mm["content"].(type) {
			case string:
				out["content"] = textBlocks(c)
			case []any:
				out["content"] = c
			default:
				add("%s: messages[%d].content must be a string or a list of blocks", where, i)
			}
			messages = append(messages, out)
		}
		if last, _ := raw[len(raw)-1].(map[string]any); last == nil || last["role"] != "user" {
			add("%s: last message must have role 'user' (assistant prefill is a 400 on Claude Mythos 5.1)", where)
		}
	}

	maxTokens := d.MaxTokens
	if v, ok := row["max_tokens"]; ok {
		n, isNum := v.(json.Number)
		i, err := int64(0), error(nil)
		if isNum {
			i, err = n.Int64()
		}
		if !isNum || err != nil || i < 1 || i > maxOutputTokens {
			add("%s: max_tokens must be an integer in 1..%d", where, maxOutputTokens)
		} else {
			maxTokens = int(i)
		}
	}

	effort := d.Effort
	if v, ok := row["effort"]; ok {
		s, isStr := v.(string)
		if !isStr || !contains(efforts, s) {
			add("%s: effort must be one of %v, got %v", where, efforts, v)
		} else {
			effort = s
		}
	}

	extra := map[string]any{}
	if v, ok := row["params"]; ok {
		m, isMap := v.(map[string]any)
		if !isMap {
			add("%s: params must be an object", where)
		} else {
			extra = m
			checkParams(extra, where, add)
		}
	}

	var system any
	switch s := row["system"].(type) {
	case nil:
		if d.System != "" {
			system = textBlocks(d.System)
		}
	case string:
		if s != "" {
			system = textBlocks(s)
		}
	case []any:
		system = s
	default:
		add("%s: system must be a string or a list of text blocks", where)
	}

	if len(problems) > 0 {
		return Request{}, problems
	}
	params := map[string]any{"model": d.Model, "max_tokens": maxTokens, "messages": messages}
	if system != nil {
		params["system"] = system
	}
	for k, v := range extra {
		params[k] = v
	}
	if effort != "" {
		oc, _ := params["output_config"].(map[string]any)
		merged := map[string]any{}
		for k, v := range oc {
			merged[k] = v
		}
		merged["effort"] = effort
		params["output_config"] = merged
	}
	return Request{CustomID: cid, Params: params}, nil
}

// loadRequests parses JSONL and reports every problem, not just the first.
func loadRequests(r io.Reader, d Defaults) ([]Request, error) {
	var out []Request
	var problems []string
	seen := map[string]int{}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 1<<20), 64<<20)
	for n := 1; sc.Scan(); n++ {
		line := bytes.TrimSpace(sc.Bytes())
		if len(line) == 0 {
			continue
		}
		where := fmt.Sprintf("line %d", n)
		dec := json.NewDecoder(bytes.NewReader(line))
		dec.UseNumber()
		var row map[string]any
		if err := dec.Decode(&row); err != nil || row == nil {
			problems = append(problems, where+": invalid JSON object")
			continue
		}
		req, p := normalizeRow(row, d, where)
		if len(p) > 0 {
			problems = append(problems, p...)
			continue
		}
		if first, dup := seen[req.CustomID]; dup {
			problems = append(problems, fmt.Sprintf("%s: duplicate custom_id %q (first on line %d)", where, req.CustomID, first))
			continue
		}
		seen[req.CustomID] = n
		out = append(out, req)
	}
	if err := sc.Err(); err != nil {
		problems = append(problems, "reading input: "+err.Error())
	}
	if len(out) == 0 && len(problems) == 0 {
		problems = append(problems, "no requests found")
	}
	if len(problems) > 0 {
		return nil, &requestError{Problems: problems}
	}
	return out, nil
}

func singlePrompt(prompt string, d Defaults) ([]Request, error) {
	req, p := normalizeRow(map[string]any{"custom_id": "prompt-1", "prompt": prompt}, d, "-prompt")
	if len(p) > 0 {
		return nil, &requestError{Problems: p}
	}
	return []Request{req}, nil
}

// chunkRequests splits into groups within the per-batch request and byte limits.
func chunkRequests(reqs []Request) ([][]Request, error) {
	var chunks [][]Request
	var cur []Request
	size := 0
	for _, r := range reqs {
		b, err := json.Marshal(r)
		if err != nil {
			return nil, err
		}
		n := len(b) + 1
		if n > maxBatchBytes {
			return nil, &requestError{Problems: []string{fmt.Sprintf("request %q alone exceeds the batch size limit", r.CustomID)}}
		}
		if len(cur) > 0 && (len(cur) >= maxRequestsPerBatch || size+n > maxBatchBytes) {
			chunks = append(chunks, cur)
			cur, size = nil, 0
		}
		cur = append(cur, r)
		size += n
	}
	if len(cur) > 0 {
		chunks = append(chunks, cur)
	}
	return chunks, nil
}
