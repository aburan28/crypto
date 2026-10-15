package main

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"
)

func load(t *testing.T, lines ...string) ([]Request, error) {
	t.Helper()
	return loadRequests(strings.NewReader(strings.Join(lines, "\n")), Defaults{Model: defaultModel, MaxTokens: defaultMaxTokens})
}

func TestPromptRowBecomesMythosRequest(t *testing.T) {
	reqs, err := load(t, `{"custom_id":"a","prompt":"hi","effort":"high","system":"be brief"}`)
	if err != nil {
		t.Fatal(err)
	}
	b, _ := json.Marshal(reqs[0].Params)
	want := `{"max_tokens":32000,"messages":[{"content":[{"text":"hi","type":"text"}],"role":"user"}],"model":"claude-mythos-5-1","output_config":{"effort":"high"},"system":[{"text":"be brief","type":"text"}]}`
	if string(b) != want {
		t.Fatalf("got  %s\nwant %s", b, want)
	}
}

func TestRejectsInvalidRows(t *testing.T) {
	cases := []struct{ row, needle string }{
		{`{"custom_id":"bad id!","prompt":"x"}`, "custom_id"},
		{`{"custom_id":"a"}`, "exactly one"},
		{`{"custom_id":"a","prompt":"x","messages":[]}`, "exactly one"},
		{`{"custom_id":"a","prompt":"x","bogus":1}`, "unknown keys"},
		{`{"custom_id":"a","prompt":"x","max_tokens":1000000}`, "max_tokens"},
		{`{"custom_id":"a","prompt":"x","max_tokens":1.5}`, "max_tokens"},
		{`{"custom_id":"a","prompt":"x","effort":"turbo"}`, "effort"},
		{`{"custom_id":"a","prompt":"x","params":{"temperature":0}}`, "temperature"},
		{`{"custom_id":"a","prompt":"x","params":{"thinking":{"type":"disabled"}}}`, "thinking"},
		{`{"custom_id":"a","prompt":"x","params":{"thinking":{"type":"enabled","budget_tokens":5000}}}`, "thinking"},
		{`{"custom_id":"a","prompt":"x","params":{"tool_choice":{"type":"any"}}}`, "tool_choice"},
		{`{"custom_id":"a","prompt":"x","params":{"model":"other"}}`, "params.model"},
		{`{"custom_id":"a","messages":[{"role":"user","content":"q"},{"role":"assistant","content":"{"}]}`, "prefill"},
		{`not json`, "invalid JSON"},
	}
	for _, c := range cases {
		_, err := load(t, c.row)
		if err == nil || !strings.Contains(err.Error()+fmt.Sprint(err.(*requestError).Problems), c.needle) {
			t.Errorf("%s: want problem containing %q, got %v", c.row, c.needle, err)
		}
	}
}

func TestReportsEveryProblem(t *testing.T) {
	_, err := load(t, `{"custom_id":"a","prompt":"x"}`, `{"custom_id":"a","prompt":"y"}`, `{nope`, ``)
	re := err.(*requestError)
	if len(re.Problems) != 2 || !strings.Contains(re.Problems[0], "duplicate") || !strings.Contains(re.Problems[1], "invalid JSON") {
		t.Fatalf("problems: %v", re.Problems)
	}
}

func TestAdaptiveThinkingAndAutoToolChoiceAllowed(t *testing.T) {
	_, err := load(t, `{"custom_id":"a","prompt":"x","params":{"thinking":{"type":"adaptive","display":"summarized"},"tool_choice":{"type":"auto"}}}`)
	if err != nil {
		t.Fatal(err)
	}
}

func TestChunking(t *testing.T) {
	reqs, _ := load(t, `{"custom_id":"a","prompt":"x"}`, `{"custom_id":"b","prompt":"x"}`, `{"custom_id":"c","prompt":"x"}`)
	old := maxRequestsPerBatch
	defer func() { maxRequestsPerBatch = old }()
	maxRequestsPerBatch = 2
	ch, _ := chunkRequests(reqs)
	if len(ch) != 2 || len(ch[0]) != 2 || len(ch[1]) != 1 {
		t.Fatalf("chunks %v", ch)
	}
	maxRequestsPerBatch = 100
	oldB := maxBatchBytes
	defer func() { maxBatchBytes = oldB }()
	maxBatchBytes = 200
	ch, _ = chunkRequests(reqs)
	if len(ch) < 2 {
		t.Fatalf("byte limit did not split: %v", ch)
	}
}

// ---- fake Batches API --------------------------------------------------

type fakeAPI struct {
	mu       sync.Mutex
	batches  map[string][]map[string]any
	creates  int
	headers  http.Header
	bodies   []map[string]any
	outcomes func(customID string, params map[string]any) string // JSON of the result object, "" = omit
	polls    int
	endAfter int // number of Get calls that report in_progress before ending
}

func (f *fakeAPI) batchJSON(id, host string) string {
	status := "ended"
	f.polls++
	if f.polls <= f.endAfter {
		status = "in_progress"
	}
	n := len(f.batches[id])
	return fmt.Sprintf(`{"id":%q,"type":"message_batch","processing_status":%q,"request_counts":{"processing":0,"succeeded":%d,"errored":0,"canceled":0,"expired":0},"created_at":"2026-10-15T00:00:00Z","expires_at":"2026-10-16T00:00:00Z","ended_at":null,"cancel_initiated_at":null,"archived_at":null,"results_url":"http://%s/v1/messages/batches/%s/results"}`, id, status, n, host, id)
}

func (f *fakeAPI) ServeHTTP(w http.ResponseWriter, r *http.Request) {
	f.mu.Lock()
	defer f.mu.Unlock()
	parts := strings.Split(strings.Trim(r.URL.Path, "/"), "/")
	w.Header().Set("content-type", "application/json")
	switch {
	case r.Method == "POST" && len(parts) == 3 && parts[2] == "batches":
		f.headers = r.Header.Clone()
		var body struct {
			Requests []map[string]any `json:"requests"`
		}
		raw, _ := io.ReadAll(r.Body)
		json.Unmarshal(raw, &body)
		f.creates++
		id := fmt.Sprintf("msgbatch_%d", f.creates)
		f.batches[id] = body.Requests
		for _, q := range body.Requests {
			f.bodies = append(f.bodies, q)
		}
		fmt.Fprint(w, f.batchJSON(id, r.Host))
	case r.Method == "GET" && len(parts) == 5 && parts[4] == "results":
		w.Header().Set("content-type", "application/binary")
		reqs := f.batches[parts[3]]
		for i := len(reqs) - 1; i >= 0; i-- { // reversed: results come back out of order
			cid := reqs[i]["custom_id"].(string)
			if res := f.outcomes(cid, reqs[i]["params"].(map[string]any)); res != "" {
				fmt.Fprintf(w, `{"custom_id":%q,"result":%s}`+"\n", cid, res)
			}
		}
	case r.Method == "GET" && len(parts) == 4:
		fmt.Fprint(w, f.batchJSON(parts[3], r.Host))
	case r.Method == "POST" && len(parts) == 5 && parts[4] == "cancel":
		j := f.batchJSON(parts[3], r.Host)
		j = strings.Replace(strings.Replace(j, `"ended"`, `"canceling"`, 1), `"in_progress"`, `"canceling"`, 1)
		fmt.Fprint(w, j)
	default:
		http.NotFound(w, r)
	}
}

func message(text, stop, extra string) string {
	return fmt.Sprintf(`{"type":"succeeded","message":{"id":"msg_1","type":"message","role":"assistant","model":"claude-mythos-5-1","content":[{"type":"thinking","thinking":"","signature":"s"},{"type":"text","text":%q}],"stop_reason":%q,"stop_sequence":null%s,"usage":{"input_tokens":10,"output_tokens":20,"cache_creation_input_tokens":0,"cache_read_input_tokens":0}}}`, text, stop, extra)
}

func startFake(t *testing.T, f *fakeAPI) {
	t.Helper()
	f.batches = map[string][]map[string]any{}
	srv := httptest.NewServer(f)
	t.Cleanup(srv.Close)
	t.Setenv("ANTHROPIC_API_KEY", "sk-ant-test")
	t.Setenv("ANTHROPIC_BASE_URL", srv.URL)
}

func runCLI(args ...string) (int, string, string) {
	var out, errb bytes.Buffer
	code := run(context.Background(), args, strings.NewReader(""), &out, &errb)
	return code, out.String(), errb.String()
}

func writeFile(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func readLines(t *testing.T, path string) []map[string]any {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var rows []map[string]any
	for _, l := range strings.Split(strings.TrimSpace(string(b)), "\n") {
		if l == "" {
			continue
		}
		var m map[string]any
		if err := json.Unmarshal([]byte(l), &m); err != nil {
			t.Fatal(err)
		}
		rows = append(rows, m)
	}
	return rows
}

func TestEndToEndRealSDK(t *testing.T) {
	f := &fakeAPI{endAfter: 1, outcomes: func(cid string, p map[string]any) string {
		return message("echo: "+cid, "end_turn", "")
	}}
	startFake(t, f)
	in := writeFile(t, "in.jsonl", `{"custom_id":"one","prompt":"alpha","effort":"high"}
{"custom_id":"two","prompt":"beta","system":"be brief","params":{"metadata":{"user_id":"u1"},"tools":[{"name":"t","description":"d","input_schema":{"type":"object","properties":{}}}]}}
`)
	out := filepath.Join(t.TempDir(), "run")
	code, stdout, stderr := runCLI("run", "-in", in, "-out", out, "-poll", "1ms")
	if code != 0 {
		t.Fatalf("exit %d\nstdout: %s\nstderr: %s", code, stdout, stderr)
	}
	// resuming must not create a second batch
	if code, so, se := runCLI("run", "-in", in, "-out", out, "-poll", "1ms"); code != 0 || !strings.Contains(so, "resuming") {
		t.Fatalf("resume exit %d: %s %s", code, so, se)
	}
	if f.creates != 1 {
		t.Fatalf("expected 1 batch create, got %d", f.creates)
	}
	if got := f.headers.Get("X-Api-Key"); got != "sk-ant-test" {
		t.Fatalf("api key header %q", got)
	}
	// the exact request bodies the SDK put on the wire
	p1 := f.bodies[0]["params"].(map[string]any)
	if p1["model"] != "claude-mythos-5-1" || p1["max_tokens"].(float64) != 32000 {
		t.Fatalf("params: %v", p1)
	}
	if oc := p1["output_config"].(map[string]any); oc["effort"] != "high" {
		t.Fatalf("output_config: %v", oc)
	}
	msg := p1["messages"].([]any)[0].(map[string]any)
	if msg["role"] != "user" || msg["content"].([]any)[0].(map[string]any)["text"] != "alpha" {
		t.Fatalf("messages: %v", msg)
	}
	p2 := f.bodies[1]["params"].(map[string]any)
	if p2["system"].([]any)[0].(map[string]any)["text"] != "be brief" || p2["metadata"].(map[string]any)["user_id"] != "u1" || len(p2["tools"].([]any)) != 1 {
		t.Fatalf("passthrough lost: %v", p2)
	}
	for _, k := range []string{"thinking", "temperature", "tool_choice", "stream"} {
		if _, bad := p1[k]; bad {
			t.Errorf("unexpected %s in body", k)
		}
	}
	rows := readLines(t, filepath.Join(out, "results.jsonl"))
	if len(rows) != 2 || rows[0]["custom_id"] != "one" || rows[0]["text"] != "echo: one" || rows[1]["custom_id"] != "two" {
		t.Fatalf("results (must be input order): %v", rows)
	}
	var s summary
	b, _ := os.ReadFile(filepath.Join(out, "summary.json"))
	json.Unmarshal(b, &s)
	if s.StatusCounts["succeeded"] != 2 || s.Tokens.OutputTokens != 40 || s.EstimatedUSD != 0.0011 {
		t.Fatalf("summary: %+v", s)
	}
	if mb, _ := os.ReadFile(filepath.Join(out, "manifest.json")); strings.Contains(string(mb), "sk-ant-test") {
		t.Fatal("manifest leaks the API key")
	}
}

func TestResultClassification(t *testing.T) {
	f := &fakeAPI{outcomes: func(cid string, p map[string]any) string {
		switch cid {
		case "ok":
			return message("answer", "end_turn", "")
		case "refused":
			return message("", "refusal", `,"stop_details":{"type":"refusal","category":"cyber","explanation":"declined"}`)
		case "cut":
			return message("partial", "max_tokens", "")
		case "bad":
			return `{"type":"errored","error":{"type":"error","error":{"type":"invalid_request_error","message":"nope"},"request_id":"req_1"}}`
		case "busy":
			return `{"type":"errored","error":{"type":"error","error":{"type":"overloaded_error","message":"busy"},"request_id":"req_2"}}`
		case "late":
			return `{"type":"expired"}`
		}
		return "" // "gone" never comes back
	}}
	startFake(t, f)
	var lines []string
	for _, id := range []string{"ok", "refused", "cut", "bad", "busy", "late", "gone"} {
		lines = append(lines, fmt.Sprintf(`{"custom_id":%q,"prompt":"p","effort":"low","params":{"metadata":{"user_id":"u"}}}`, id))
	}
	in := writeFile(t, "in.jsonl", strings.Join(lines, "\n"))
	out := filepath.Join(t.TempDir(), "run")
	code, _, stderr := runCLI("run", "-in", in, "-out", out, "-poll", "1ms")
	if code != exitFailures {
		t.Fatalf("exit %d, stderr %s", code, stderr)
	}
	status := map[string]string{}
	rows := readLines(t, filepath.Join(out, "results.jsonl"))
	for _, r := range rows {
		status[r["custom_id"].(string)] = r["status"].(string)
	}
	want := map[string]string{"ok": "succeeded", "refused": "refused", "cut": "succeeded", "bad": "errored", "busy": "errored", "late": "expired", "gone": "missing"}
	for k, v := range want {
		if status[k] != v {
			t.Errorf("%s: status %q, want %q", k, status[k], v)
		}
	}
	if sd := rows[1]["stop_details"].(map[string]any); sd["category"] != "cyber" {
		t.Errorf("stop_details %v", sd)
	}
	if rows[2]["truncated"] != true {
		t.Errorf("truncation not flagged: %v", rows[2])
	}
	if e := rows[3]["error"].(map[string]any); e["type"] != "invalid_request_error" || e["request_id"] != "req_1" {
		t.Errorf("error %v", e)
	}
	// retryable: busy, late, gone — not the invalid request, not the refusal
	retry := readLines(t, filepath.Join(out, "retryable.jsonl"))
	var ids []string
	for _, r := range retry {
		ids = append(ids, r["custom_id"].(string))
	}
	if strings.Join(ids, ",") != "busy,late,gone" {
		t.Fatalf("retryable %v", ids)
	}
	// and the retry file is itself valid input that preserves effort and params
	f2, _ := os.Open(filepath.Join(out, "retryable.jsonl"))
	defer f2.Close()
	again, err := loadRequests(f2, Defaults{Model: defaultModel, MaxTokens: defaultMaxTokens})
	if err != nil || len(again) != 3 {
		t.Fatalf("retry file not reloadable: %v", err)
	}
	if again[0].Params["output_config"].(map[string]any)["effort"] != "low" || again[0].Params["metadata"] == nil {
		t.Fatalf("retry lost fields: %v", again[0].Params)
	}
}

func TestSubmitRefusesSecondSubmissionAndMultiBatch(t *testing.T) {
	f := &fakeAPI{outcomes: func(string, map[string]any) string { return "" }}
	startFake(t, f)
	old := maxRequestsPerBatch
	defer func() { maxRequestsPerBatch = old }()
	maxRequestsPerBatch = 2
	in := writeFile(t, "in.jsonl", `{"custom_id":"a","prompt":"x"}
{"custom_id":"b","prompt":"x"}
{"custom_id":"c","prompt":"x"}`)
	out := filepath.Join(t.TempDir(), "run")
	if code, _, se := runCLI("submit", "-in", in, "-out", out); code != 0 {
		t.Fatalf("submit: %d %s", code, se)
	}
	if f.creates != 2 {
		t.Fatalf("expected 2 batches, got %d", f.creates)
	}
	code, _, se := runCLI("submit", "-in", in, "-out", out)
	if code != exitUsage || !strings.Contains(se, "already has a submission") || f.creates != 2 {
		t.Fatalf("second submit: code %d creates %d %s", code, f.creates, se)
	}
	// run with different requests into the same directory is refused too
	other := writeFile(t, "other.jsonl", `{"custom_id":"z","prompt":"different"}`)
	if code, _, se := runCLI("run", "-in", other, "-out", out, "-poll", "1ms"); code != exitUsage || !strings.Contains(se, "different requests") {
		t.Fatalf("run with other input: %d %s", code, se)
	}
}

func TestCollectBeforeEnd(t *testing.T) {
	f := &fakeAPI{endAfter: 1000, outcomes: func(string, map[string]any) string { return "" }}
	startFake(t, f)
	in := writeFile(t, "in.jsonl", `{"custom_id":"a","prompt":"x"}`)
	out := filepath.Join(t.TempDir(), "run")
	runCLI("submit", "-in", in, "-out", out)
	if code, _, se := runCLI("collect", "-out", out); code != exitUsage || !strings.Contains(se, "not ended") {
		t.Fatalf("collect early: %d %s", code, se)
	}
	if code, so, _ := runCLI("status", "-out", out); code != 0 || !strings.Contains(so, "in_progress") {
		t.Fatalf("status: %d %s", code, so)
	}
	if code, so, _ := runCLI("cancel", "-out", out); code != 0 || !strings.Contains(so, "canceling") {
		t.Fatalf("cancel: %d %s", code, so)
	}
}

func TestWaitTimeoutAndInterrupt(t *testing.T) {
	f := &fakeAPI{endAfter: 1000, outcomes: func(string, map[string]any) string { return "" }}
	startFake(t, f)
	in := writeFile(t, "in.jsonl", `{"custom_id":"a","prompt":"x"}`)
	out := filepath.Join(t.TempDir(), "run")
	runCLI("submit", "-in", in, "-out", out)
	if code, _, se := runCLI("wait", "-out", out, "-poll", "1ms", "-timeout", "20ms"); code != exitUsage || !strings.Contains(se, "keep running server-side") {
		t.Fatalf("timeout: %d %s", code, se)
	}
	ctx, cancel := context.WithCancel(context.Background())
	defer cancel()
	time.AfterFunc(30*time.Millisecond, cancel) // what Ctrl-C does via signal.NotifyContext
	var so, se bytes.Buffer
	if code := run(ctx, []string{"wait", "-out", out, "-poll", "1h"}, strings.NewReader(""), &so, &se); code != exitInterrupted {
		t.Fatalf("interrupt: %d %s", code, se.String())
	}
}

func TestMissingKeyAndAPIErrorHints(t *testing.T) {
	t.Setenv("ANTHROPIC_API_KEY", "")
	t.Setenv("ANTHROPIC_AUTH_TOKEN", "")
	if code, _, se := runCLI("submit", "-prompt", "hi", "-out", t.TempDir()+"/o"); code != exitUsage || !strings.Contains(se, "ANTHROPIC_API_KEY") {
		t.Fatalf("missing key: %d %s", code, se)
	}
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		w.WriteHeader(400)
		fmt.Fprint(w, `{"type":"error","error":{"type":"invalid_request_error","message":"In order to access this model, your organization or workspace must have data retention enabled."}}`)
	}))
	defer srv.Close()
	t.Setenv("ANTHROPIC_API_KEY", "sk-ant-test")
	t.Setenv("ANTHROPIC_BASE_URL", srv.URL)
	code, _, se := runCLI("submit", "-prompt", "hi", "-out", t.TempDir()+"/o")
	if code != exitAPI || !strings.Contains(se, "30-day data retention") {
		t.Fatalf("retention hint: %d %s", code, se)
	}
}

func TestValidateNeedsNoKey(t *testing.T) {
	t.Setenv("ANTHROPIC_API_KEY", "")
	in := writeFile(t, "in.jsonl", `{"custom_id":"a","prompt":"x"}`)
	if code, so, _ := runCLI("validate", "-in", in); code != 0 || !strings.Contains(so, "claude-mythos-5-1") {
		t.Fatalf("validate: %d %s", code, so)
	}
	bad := writeFile(t, "bad.jsonl", `{"custom_id":"a","prompt":"x","params":{"temperature":1}}`)
	if code, _, _ := runCLI("validate", "-in", bad); code != exitUsage {
		t.Fatalf("validate bad: %d", code)
	}
	if code, _, se := runCLI("validate"); code != exitUsage || !strings.Contains(se, "exactly one of -in or -prompt") {
		t.Fatalf("no input: %d %s", code, se)
	}
}
