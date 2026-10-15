package main

// Submit, wait for and collect Message Batches. State lives in an output
// directory:
//
//	manifest.json    model, hash of the requests, batch ids (never credentials)
//	requests.jsonl   the normalized requests as submitted
//	results.jsonl    one row per request, in input order (after collect)
//	retryable.jsonl  input rows worth resubmitting
//	summary.json     counts, token usage, estimated cost
//
// submit refuses a directory that already holds a submission, so a re-run
// cannot bill a second batch; run resumes from the manifest instead.

import (
	"context"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"github.com/anthropics/anthropic-sdk-go"
	"github.com/anthropics/anthropic-sdk-go/option"
)

// Batch list price per MTok for Claude Mythos 5.1 (50% of $10 / $50). Estimate
// only: cache reads and writes are not priced.
const (
	batchInputPerMTok  = 5.00
	batchOutputPerMTok = 25.00
	manifestName       = "manifest.json"
	manifestSchema     = 1
)

var errState = errors.New("state")

func stateErr(format string, a ...any) error {
	return fmt.Errorf("%w: %s", errState, fmt.Sprintf(format, a...))
}

type batchEntry struct {
	ID            string `json:"id"`
	NRequests     int    `json:"n_requests"`
	FirstCustomID string `json:"first_custom_id"`
}

type manifest struct {
	Schema      int          `json:"schema"`
	Model       string       `json:"model"`
	CreatedAt   string       `json:"created_at"`
	InputSHA256 string       `json:"input_sha256"`
	NRequests   int          `json:"n_requests"`
	Batches     []batchEntry `json:"batches"`
}

type logf func(format string, a ...any)

func newClient() (anthropic.Client, error) {
	if os.Getenv("ANTHROPIC_API_KEY") == "" && os.Getenv("ANTHROPIC_AUTH_TOKEN") == "" {
		return anthropic.Client{}, errors.New("no Anthropic credentials found: export ANTHROPIC_API_KEY")
	}
	// ANTHROPIC_BASE_URL is honoured by the SDK for proxies.
	return anthropic.NewClient(option.WithMaxRetries(5)), nil
}

func now() string { return time.Now().UTC().Format(time.RFC3339) }

func jsonl[T any](rows []T) ([]byte, error) {
	var b []byte
	for _, r := range rows {
		line, err := json.Marshal(r)
		if err != nil {
			return nil, err
		}
		b = append(append(b, line...), '\n')
	}
	return b, nil
}

func atomicWrite(path string, data []byte) error {
	tmp := path + ".tmp"
	if err := os.WriteFile(tmp, data, 0o644); err != nil {
		return err
	}
	return os.Rename(tmp, path)
}

func writeJSON(path string, v any) error {
	b, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return err
	}
	return atomicWrite(path, append(b, '\n'))
}

func requestDigest(reqs []Request) (string, error) {
	b, err := jsonl(reqs)
	if err != nil {
		return "", err
	}
	sum := sha256.Sum256(b)
	return hex.EncodeToString(sum[:]), nil
}

func loadManifest(out string) (manifest, error) {
	var m manifest
	b, err := os.ReadFile(filepath.Join(out, manifestName))
	if err != nil {
		return m, stateErr("%s not found; run submit first", filepath.Join(out, manifestName))
	}
	return m, json.Unmarshal(b, &m)
}

// toBatchRequest turns a stored request into the SDK type. The SDK param
// structs unmarshal from the same JSON the API takes, so tools, metadata and
// other per-row params pass straight through.
func toBatchRequest(r Request) (anthropic.MessageBatchNewParamsRequest, error) {
	raw, err := json.Marshal(r.Params)
	if err != nil {
		return anthropic.MessageBatchNewParamsRequest{}, err
	}
	var p anthropic.MessageBatchNewParamsRequestParams
	if err := json.Unmarshal(raw, &p); err != nil {
		return anthropic.MessageBatchNewParamsRequest{}, fmt.Errorf("request %q: %w", r.CustomID, err)
	}
	return anthropic.MessageBatchNewParamsRequest{CustomID: r.CustomID, Params: p}, nil
}

func submit(ctx context.Context, c anthropic.Client, reqs []Request, out, model string, log logf) (manifest, error) {
	if err := os.MkdirAll(out, 0o755); err != nil {
		return manifest{}, err
	}
	if _, err := os.Stat(filepath.Join(out, manifestName)); err == nil {
		return manifest{}, stateErr("%s exists: this directory already has a submission (use wait/collect, or pick a new -out)", filepath.Join(out, manifestName))
	}
	chunks, err := chunkRequests(reqs)
	if err != nil {
		return manifest{}, err
	}
	// Convert everything first: a bad row must fail before any batch is billed.
	typed := make([][]anthropic.MessageBatchNewParamsRequest, len(chunks))
	for i, ch := range chunks {
		for _, r := range ch {
			t, err := toBatchRequest(r)
			if err != nil {
				return manifest{}, &requestError{Problems: []string{err.Error()}}
			}
			typed[i] = append(typed[i], t)
		}
	}
	digest, err := requestDigest(reqs)
	if err != nil {
		return manifest{}, err
	}
	m := manifest{Schema: manifestSchema, Model: model, CreatedAt: now(), InputSHA256: digest, NRequests: len(reqs)}
	body, err := jsonl(reqs)
	if err != nil {
		return m, err
	}
	if err := atomicWrite(filepath.Join(out, "requests.jsonl"), body); err != nil {
		return m, err
	}
	for i, ch := range chunks {
		b, err := c.Messages.Batches.New(ctx, anthropic.MessageBatchNewParams{Requests: typed[i]})
		if err != nil {
			return m, err
		}
		m.Batches = append(m.Batches, batchEntry{ID: b.ID, NRequests: len(ch), FirstCustomID: ch[0].CustomID})
		// Persist after every batch: an interrupted multi-batch submit stays recoverable.
		if err := writeJSON(filepath.Join(out, manifestName), m); err != nil {
			return m, err
		}
		log("submitted batch %d/%d: %s (%d requests)", i+1, len(chunks), b.ID, len(ch))
	}
	return m, nil
}

type statusRow struct {
	ID               string `json:"id"`
	ProcessingStatus string `json:"processing_status"`
	Processing       int64  `json:"processing"`
	Succeeded        int64  `json:"succeeded"`
	Errored          int64  `json:"errored"`
	Canceled         int64  `json:"canceled"`
	Expired          int64  `json:"expired"`
}

func status(ctx context.Context, c anthropic.Client, out string) ([]statusRow, error) {
	m, err := loadManifest(out)
	if err != nil {
		return nil, err
	}
	var rows []statusRow
	for _, e := range m.Batches {
		b, err := c.Messages.Batches.Get(ctx, e.ID, anthropic.MessageBatchGetParams{})
		if err != nil {
			return nil, err
		}
		n := b.RequestCounts
		rows = append(rows, statusRow{b.ID, string(b.ProcessingStatus), n.Processing, n.Succeeded, n.Errored, n.Canceled, n.Expired})
	}
	return rows, nil
}

var errTimeout = errors.New("timeout")

// wait is the loop: poll every batch until all have ended. Batches finish in
// 24 h at most; they keep running server-side if this process stops.
func wait(ctx context.Context, c anthropic.Client, out string, poll, timeout time.Duration, log logf) error {
	start := time.Now()
	for {
		rows, err := status(ctx, c, out)
		if err != nil {
			return err
		}
		done := 0
		var parts []string
		for _, r := range rows {
			if r.ProcessingStatus == "ended" {
				done++
			}
			parts = append(parts, fmt.Sprintf("%s: %s (ok %d, err %d, pending %d)", r.ID, r.ProcessingStatus, r.Succeeded, r.Errored, r.Processing))
		}
		log("[%s] %d/%d batches ended; %s", now(), done, len(rows), strings.Join(parts, ", "))
		if done == len(rows) {
			return nil
		}
		if time.Since(start) > timeout {
			return fmt.Errorf("%w: batches still running after %.1f h; they keep running server-side - resume with: wait -out %s", errTimeout, timeout.Hours(), out)
		}
		select {
		case <-ctx.Done():
			return ctx.Err()
		case <-time.After(poll):
		}
	}
}

type usage struct {
	InputTokens         int64 `json:"input_tokens"`
	OutputTokens        int64 `json:"output_tokens"`
	CacheCreationTokens int64 `json:"cache_creation_input_tokens"`
	CacheReadTokens     int64 `json:"cache_read_input_tokens"`
}

type errInfo struct {
	Type      string `json:"type"`
	Message   string `json:"message"`
	RequestID string `json:"request_id,omitempty"`
}

type stopDetails struct {
	Category    string `json:"category"`
	Explanation string `json:"explanation"`
}

type resultRow struct {
	CustomID    string       `json:"custom_id"`
	Status      string       `json:"status"`
	StopReason  string       `json:"stop_reason,omitempty"`
	Model       string       `json:"model,omitempty"`
	Text        *string      `json:"text,omitempty"`
	Usage       *usage       `json:"usage,omitempty"`
	Truncated   bool         `json:"truncated,omitempty"`
	StopDetails *stopDetails `json:"stop_details,omitempty"`
	Error       *errInfo     `json:"error,omitempty"`
}

func toResultRow(item anthropic.MessageBatchIndividualResponse) resultRow {
	row := resultRow{CustomID: item.CustomID}
	switch v := item.Result.AsAny().(type) {
	case anthropic.MessageBatchSucceededResult:
		msg := v.Message
		var text strings.Builder
		for _, b := range msg.Content {
			if tb, ok := b.AsAny().(anthropic.TextBlock); ok {
				text.WriteString(tb.Text)
			}
		}
		t := text.String()
		row.Status, row.StopReason, row.Model, row.Text = "succeeded", string(msg.StopReason), string(msg.Model), &t
		row.Usage = &usage{msg.Usage.InputTokens, msg.Usage.OutputTokens, msg.Usage.CacheCreationInputTokens, msg.Usage.CacheReadInputTokens}
		switch msg.StopReason {
		case anthropic.StopReasonRefusal:
			// Safeguards declined: HTTP 200, so the API calls it "succeeded".
			row.Status = "refused"
			row.StopDetails = &stopDetails{string(msg.StopDetails.Category), msg.StopDetails.Explanation}
		case anthropic.StopReasonMaxTokens:
			row.Truncated = true
		}
	case anthropic.MessageBatchErroredResult:
		row.Status = "errored"
		row.Error = &errInfo{Type: v.Error.Error.Type, Message: v.Error.Error.Message, RequestID: v.Error.RequestID}
	case anthropic.MessageBatchCanceledResult:
		row.Status = "canceled"
	case anthropic.MessageBatchExpiredResult:
		row.Status = "expired"
	default:
		row.Status = "unknown"
	}
	return row
}

type summary struct {
	Model        string         `json:"model"`
	CollectedAt  string         `json:"collected_at"`
	Batches      []string       `json:"batches"`
	NRequests    int            `json:"n_requests"`
	StatusCounts map[string]int `json:"status_counts"`
	Truncated    int            `json:"truncated_max_tokens"`
	Retryable    int            `json:"retryable"`
	Tokens       usage          `json:"tokens"`
	EstimatedUSD float64        `json:"estimated_cost_usd_input_output_only"`
	CostNote     string         `json:"cost_note"`
}

func readRequests(out string) ([]Request, error) {
	f, err := os.Open(filepath.Join(out, "requests.jsonl"))
	if err != nil {
		return nil, stateErr("requests.jsonl missing in %s", out)
	}
	defer f.Close()
	// requests.jsonl holds normalized {custom_id, params}; decode it directly.
	var reqs []Request
	dec := json.NewDecoder(f)
	for dec.More() {
		var r Request
		if err := dec.Decode(&r); err != nil {
			return nil, err
		}
		reqs = append(reqs, r)
	}
	return reqs, nil
}

func collect(ctx context.Context, c anthropic.Client, out string, log logf) (summary, error) {
	m, err := loadManifest(out)
	if err != nil {
		return summary{}, err
	}
	reqs, err := readRequests(out)
	if err != nil {
		return summary{}, err
	}
	byID := map[string]resultRow{}
	var pending []string
	for _, e := range m.Batches {
		b, err := c.Messages.Batches.Get(ctx, e.ID, anthropic.MessageBatchGetParams{})
		if err != nil {
			return summary{}, err
		}
		if b.ProcessingStatus != anthropic.MessageBatchProcessingStatusEnded {
			pending = append(pending, e.ID)
			continue
		}
		// Results arrive in any order: key by custom_id, never by position.
		stream := c.Messages.Batches.ResultsStreaming(ctx, e.ID, anthropic.MessageBatchResultsParams{})
		for stream.Next() {
			row := toResultRow(stream.Current())
			byID[row.CustomID] = row
		}
		err = stream.Err()
		stream.Close()
		if err != nil {
			return summary{}, err
		}
	}
	if len(pending) > 0 {
		return summary{}, stateErr("batches not ended yet: %s; run wait first", strings.Join(pending, ", "))
	}

	s := summary{Model: m.Model, CollectedAt: now(), NRequests: len(reqs), StatusCounts: map[string]int{}}
	rows := make([]resultRow, 0, len(reqs))
	var retry []map[string]any
	for _, r := range reqs {
		row, ok := byID[r.CustomID]
		if !ok {
			row = resultRow{CustomID: r.CustomID, Status: "missing"}
		}
		rows = append(rows, row)
		s.StatusCounts[row.Status]++
		if row.Truncated {
			s.Truncated++
		}
		if row.Usage != nil {
			s.Tokens.InputTokens += row.Usage.InputTokens
			s.Tokens.OutputTokens += row.Usage.OutputTokens
			s.Tokens.CacheCreationTokens += row.Usage.CacheCreationTokens
			s.Tokens.CacheReadTokens += row.Usage.CacheReadTokens
		}
		// invalid_request_error rows fail again unchanged; everything else may be retried.
		switch {
		case row.Status == "expired", row.Status == "canceled", row.Status == "missing",
			row.Status == "errored" && row.Error != nil && row.Error.Type != "invalid_request_error":
			retry = append(retry, retryRow(r))
		}
	}
	s.Retryable = len(retry)
	s.Batches = nil
	for _, e := range m.Batches {
		s.Batches = append(s.Batches, e.ID)
	}
	s.EstimatedUSD = float64(int((float64(s.Tokens.InputTokens)*batchInputPerMTok+float64(s.Tokens.OutputTokens)*batchOutputPerMTok)/1e6*10000+0.5)) / 10000
	s.CostNote = fmt.Sprintf("batch list price $%.0f/$%.0f per MTok; cache tokens not priced; confirm against the Console usage page", batchInputPerMTok, batchOutputPerMTok)

	for name, v := range map[string]any{"results.jsonl": rows, "retryable.jsonl": retry} {
		var b []byte
		switch x := v.(type) {
		case []resultRow:
			b, err = jsonl(x)
		case []map[string]any:
			b, err = jsonl(x)
		}
		if err != nil {
			return s, err
		}
		if err := atomicWrite(filepath.Join(out, name), b); err != nil {
			return s, err
		}
	}
	if err := writeJSON(filepath.Join(out, "summary.json"), s); err != nil {
		return s, err
	}
	log("collected %d results: %s; estimated $%.2f", len(rows), countsString(s.StatusCounts), s.EstimatedUSD)
	return s, nil
}

func countsString(c map[string]int) string {
	keys := make([]string, 0, len(c))
	for k := range c {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	parts := make([]string, 0, len(keys))
	for _, k := range keys {
		parts = append(parts, fmt.Sprintf("%s=%d", k, c[k]))
	}
	return strings.Join(parts, " ")
}

// retryRow is the inverse of normalizeRow: a valid input row for resubmission
// (the model comes from the next run's flags).
func retryRow(r Request) map[string]any {
	row := map[string]any{"custom_id": r.CustomID, "messages": r.Params["messages"], "max_tokens": r.Params["max_tokens"]}
	if s, ok := r.Params["system"]; ok {
		row["system"] = s
	}
	extra := map[string]any{}
	for k, v := range r.Params {
		switch k {
		case "model", "messages", "max_tokens", "system", "output_config":
		default:
			extra[k] = v
		}
	}
	if oc, ok := r.Params["output_config"].(map[string]any); ok {
		rest := map[string]any{}
		for k, v := range oc {
			if k == "effort" {
				row["effort"] = v
			} else {
				rest[k] = v
			}
		}
		if len(rest) > 0 {
			extra["output_config"] = rest
		}
	}
	if len(extra) > 0 {
		row["params"] = extra
	}
	return row
}

func cancelAll(ctx context.Context, c anthropic.Client, out string, log logf) error {
	m, err := loadManifest(out)
	if err != nil {
		return err
	}
	for _, e := range m.Batches {
		b, err := c.Messages.Batches.Cancel(ctx, e.ID, anthropic.MessageBatchCancelParams{})
		if err != nil {
			return err
		}
		log("%s: %s", e.ID, b.ProcessingStatus)
	}
	return nil
}

// describeAPIError gives actionable text for the failures this tool is likely to meet.
func describeAPIError(err error) (string, bool) {
	var apiErr *anthropic.Error
	if !errors.As(err, &apiErr) {
		return "", false
	}
	hint := map[int]string{
		401: "ANTHROPIC_API_KEY is missing, malformed or revoked.",
		403: "this key may lack access to claude-mythos-5-1 (Project Glasswing participants only) or to the Batches API.",
		404: "model or batch not found; claude-mythos-5-1 is hidden from keys outside Project Glasswing, and batch results expire after 29 days.",
		413: "batch exceeds 256 MB; split the input file.",
		429: "rate limited; the Mythos models have their own rate-limit pool. Retry later or submit fewer requests.",
	}[apiErr.StatusCode]
	if apiErr.StatusCode == 400 {
		if strings.Contains(strings.ToLower(err.Error()), "data retention") {
			hint = "Mythos requires 30-day data retention on the organization/workspace (not available under zero data retention)."
		} else {
			hint = "the API rejected the request body; fix the input rows."
		}
	}
	msg := fmt.Sprintf("API error (HTTP %d): %s", apiErr.StatusCode, err)
	if hint != "" {
		msg += "\n  hint: " + hint
	}
	return msg, true
}
