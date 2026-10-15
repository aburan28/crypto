// mythos-batch runs prompts through the Anthropic Message Batches API with
// Claude Mythos 5.1: submit, poll until the batch ends, write the results.
//
//	export ANTHROPIC_API_KEY=sk-ant-...
//	mythos-batch run -in requests.jsonl -out runs/my-run
//
// See README.md for the input format and the output directory.
package main

import (
	"context"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"os/signal"
	"strings"
	"syscall"
	"time"

	"github.com/anthropics/anthropic-sdk-go"
	"github.com/anthropics/anthropic-sdk-go/packages/param"
)

const (
	exitOK          = 0
	exitUsage       = 2
	exitAPI         = 3
	exitFailures    = 4
	exitInterrupted = 130
)

const usage_ = `usage: mythos-batch <command> [flags]

commands:
  validate  check the input locally (no key, no cost)
  run       submit (or resume) + wait + collect
  submit    create the batch(es) and return
  status    show batch progress
  wait      block until the batches end
  collect   download results
  cancel    cancel the batches
  list      recent batches on this workspace

input flags (validate, submit, run):
  -in FILE      JSONL, one request per line ("-" for stdin)
  -prompt TEXT  a single prompt instead (custom_id prompt-1), for smoke tests
  -model ID     default $MYTHOS_MODEL or ` + defaultModel + `
  -max-tokens N default for rows without max_tokens (thinking counts)
  -effort LEVEL default output_config.effort (low|medium|high|xhigh|max)
  -system TEXT  default system prompt for rows without one
other flags: -out DIR (all but validate, list); -poll DURATION, -timeout DURATION (run, wait)

environment: ANTHROPIC_API_KEY (required), ANTHROPIC_BASE_URL, MYTHOS_MODEL, MYTHOS_MAX_TOKENS
exit codes: 0 ok, 2 bad input/usage/credentials, 3 API error, 4 finished with non-succeeded rows, 130 interrupted
`

type options struct {
	in, prompt, out string
	defaults        Defaults
	poll, timeout   time.Duration
	limit           int
}

func parseFlags(cmd string, args []string, stderr io.Writer) (options, error) {
	d, err := defaultsFromEnv()
	if err != nil {
		return options{}, err
	}
	o := options{defaults: d}
	fs := flag.NewFlagSet(cmd, flag.ContinueOnError)
	fs.SetOutput(stderr)
	needsInput := cmd == "validate" || cmd == "submit" || cmd == "run"
	if needsInput {
		fs.StringVar(&o.in, "in", "", "")
		fs.StringVar(&o.prompt, "prompt", "", "")
		fs.StringVar(&o.defaults.Model, "model", d.Model, "")
		fs.IntVar(&o.defaults.MaxTokens, "max-tokens", d.MaxTokens, "")
		fs.StringVar(&o.defaults.Effort, "effort", "", "")
		fs.StringVar(&o.defaults.System, "system", "", "")
	}
	if cmd != "validate" && cmd != "list" {
		fs.StringVar(&o.out, "out", "", "")
	}
	if cmd == "run" || cmd == "wait" {
		fs.DurationVar(&o.poll, "poll", time.Minute, "")
		fs.DurationVar(&o.timeout, "timeout", 25*time.Hour, "")
	}
	if cmd == "list" {
		fs.IntVar(&o.limit, "limit", 20, "")
	}
	if err := fs.Parse(args); err != nil {
		return o, err
	}
	if fs.NArg() > 0 {
		return o, fmt.Errorf("unexpected argument %q", fs.Arg(0))
	}
	if needsInput {
		if (o.in == "") == (o.prompt == "") {
			return o, errors.New("give exactly one of -in or -prompt")
		}
		if p := o.defaults.check(); len(p) > 0 {
			return o, errors.New(strings.Join(p, "; "))
		}
	}
	if cmd != "validate" && cmd != "list" && o.out == "" {
		return o, errors.New("-out DIR is required")
	}
	return o, nil
}

func loadInput(o options, stdin io.Reader) ([]Request, error) {
	if o.prompt != "" {
		return singlePrompt(o.prompt, o.defaults)
	}
	if o.in == "-" {
		return loadRequests(stdin, o.defaults)
	}
	f, err := os.Open(o.in)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return loadRequests(f, o.defaults)
}

func run(ctx context.Context, args []string, stdin io.Reader, stdout, stderr io.Writer) int {
	if len(args) == 0 || args[0] == "-h" || args[0] == "-help" || args[0] == "help" {
		fmt.Fprint(stderr, usage_)
		return exitUsage
	}
	cmd := args[0]
	switch cmd {
	case "validate", "run", "submit", "status", "wait", "collect", "cancel", "list":
	default:
		fmt.Fprintf(stderr, "unknown command %q\n\n%s", cmd, usage_)
		return exitUsage
	}
	o, err := parseFlags(cmd, args[1:], stderr)
	if err != nil {
		if !errors.Is(err, flag.ErrHelp) {
			fmt.Fprintln(stderr, "error:", err)
		}
		return exitUsage
	}
	log := func(format string, a ...any) { fmt.Fprintf(stdout, format+"\n", a...) }

	fail := func(err error) int {
		var re *requestError
		switch {
		case errors.As(err, &re):
			fmt.Fprintln(stderr, "invalid input:\n  "+strings.Join(re.Problems, "\n  "))
			return exitUsage
		case errors.Is(err, errState), errors.Is(err, errTimeout), errors.Is(err, os.ErrNotExist):
			fmt.Fprintln(stderr, "error:", err)
			return exitUsage
		case errors.Is(err, context.Canceled):
			fmt.Fprintln(stderr, "\ninterrupted; submitted batches keep running server-side - resume with wait/collect (or re-run run) on the same -out")
			return exitInterrupted
		}
		if msg, ok := describeAPIError(err); ok {
			fmt.Fprintln(stderr, msg)
			return exitAPI
		}
		fmt.Fprintln(stderr, "error:", err)
		return exitAPI
	}

	var reqs []Request
	if cmd == "validate" || cmd == "submit" || cmd == "run" {
		if reqs, err = loadInput(o, stdin); err != nil {
			return fail(err)
		}
	}
	if cmd == "validate" {
		chunks, err := chunkRequests(reqs)
		if err != nil {
			return fail(err)
		}
		log("ok: %d requests, model %s, %d batch(es)", len(reqs), o.defaults.Model, len(chunks))
		return exitOK
	}

	client, err := newClient()
	if err != nil {
		fmt.Fprintln(stderr, "error:", err)
		return exitUsage
	}

	switch cmd {
	case "submit", "run":
		_, statErr := os.Stat(o.out + "/" + manifestName)
		if cmd == "run" && statErr == nil {
			m, err := loadManifest(o.out)
			if err != nil {
				return fail(err)
			}
			digest, err := requestDigest(reqs)
			if err != nil {
				return fail(err)
			}
			if m.InputSHA256 != digest {
				return fail(stateErr("%s holds a submission for different requests; pick a new -out", o.out))
			}
			log("resuming %s (already submitted)", o.out)
		} else if _, err := submit(ctx, client, reqs, o.out, o.defaults.Model, log); err != nil {
			return fail(err)
		}
		if cmd == "submit" {
			log("submitted; continue with: mythos-batch wait -out %s && mythos-batch collect -out %s", o.out, o.out)
			return exitOK
		}
		fallthrough
	case "wait":
		if err := wait(ctx, client, o.out, o.poll, o.timeout, log); err != nil {
			return fail(err)
		}
		if cmd == "wait" {
			return exitOK
		}
		fallthrough
	case "collect":
		s, err := collect(ctx, client, o.out, log)
		if err != nil {
			return fail(err)
		}
		return resultsExit(s, stderr)
	case "status":
		rows, err := status(ctx, client, o.out)
		if err != nil {
			return fail(err)
		}
		for _, r := range rows {
			fmt.Fprintf(stdout, "%s %s processing=%d succeeded=%d errored=%d canceled=%d expired=%d\n",
				r.ID, r.ProcessingStatus, r.Processing, r.Succeeded, r.Errored, r.Canceled, r.Expired)
		}
	case "cancel":
		if err := cancelAll(ctx, client, o.out, log); err != nil {
			return fail(err)
		}
	case "list":
		it := client.Messages.Batches.ListAutoPaging(ctx, anthropic.MessageBatchListParams{Limit: param.NewOpt(int64(min(o.limit, 100)))})
		for n := 0; it.Next() && n < o.limit; n++ {
			b := it.Current()
			fmt.Fprintln(stdout, b.ID, b.ProcessingStatus, b.CreatedAt.Format(time.RFC3339))
		}
		if err := it.Err(); err != nil {
			return fail(err)
		}
	}
	return exitOK
}

func resultsExit(s summary, stderr io.Writer) int {
	bad := map[string]int{}
	for k, v := range s.StatusCounts {
		if k != "succeeded" {
			bad[k] = v
		}
	}
	if len(bad) > 0 {
		fmt.Fprintf(stderr, "not all requests succeeded: %s (see results.jsonl; retryable.jsonl lists resubmittable rows)\n", countsString(bad))
		return exitFailures
	}
	if s.Truncated > 0 {
		fmt.Fprintf(stderr, "warning: %d responses hit max_tokens (thinking counts toward it); raise -max-tokens and resubmit those\n", s.Truncated)
	}
	return exitOK
}

func main() {
	ctx, stop := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer stop()
	os.Exit(run(ctx, os.Args[1:], os.Stdin, os.Stdout, os.Stderr))
}
