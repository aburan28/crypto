package main

import (
	"bytes"
	"context"
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"time"

	modal "github.com/modal-labs/modal-client/go"
)

type fileEntry struct {
	Local  string `json:"local"`
	Remote string `json:"remote"`
	SHA256 string `json:"sha256"`
	Bytes  int64  `json:"bytes"`
}

type plan struct {
	Kind                    string            `json:"kind"`
	SDKVersion              string            `json:"modal_go_sdk_version"`
	GoVersion               string            `json:"go_version"`
	GoModSHA256             string            `json:"go_mod_sha256"`
	GoSumSHA256             string            `json:"go_sum_sha256"`
	LaunchSourceSHA256      string            `json:"launch_source_sha256"`
	AnalyzerSourceSHA256    string            `json:"analyzer_source_sha256"`
	GuestSourceSHA256       map[string]string `json:"guest_source_sha256"`
	GuestBuildReceiptSHA256 string            `json:"guest_build_receipt_sha256"`
	PreviousReceiptSHA256   string            `json:"previous_receipt_sha256"`
	BundleSHA256            string            `json:"bundle_sha256"`
	BinarySHA256            map[string]string `json:"binary_sha256"`
	GuestSHA256             string            `json:"guest_sha256"`
	Files                   []fileEntry       `json:"files"`
	TransferredBytes        int64             `json:"transferred_bytes"`
	SandboxCreationAttempts int               `json:"sandbox_creation_attempts"`
	CPUPhysicalCores        int               `json:"cpu_physical_cores_request_and_limit"`
	MemoryMiB               int               `json:"memory_mib_request_and_limit"`
	LogicalExecutionCPUsMax int               `json:"logical_execution_cpus_max"`
	SandboxTimeoutSeconds   int               `json:"sandbox_timeout_seconds"`
	WorkloadDeadlineSeconds int               `json:"workload_deadline_seconds"`
	PerRunTimeoutSeconds    int               `json:"per_run_timeout_seconds"`
	MinimumRunSeconds       int               `json:"minimum_run_seconds"`
	BlockNetwork            bool              `json:"block_network"`
	VMRuntime               bool              `json:"vm_runtime"`
	Region                  *string           `json:"region"`
	Image                   string            `json:"image"`
	AppName                 string            `json:"app_name"`
	PublishedBaseComputeUSD string            `json:"published_base_compute_usd_for_timeout"`
}

type attempt struct {
	Kind               string      `json:"kind"`
	FreezeSHA256       string      `json:"freeze_sha256"`
	StartedUnixNS      int64       `json:"started_unix_ns"`
	SandboxCreateCalls int         `json:"sandbox_create_calls"`
	CPUPhysicalCores   int         `json:"cpu_physical_cores"`
	MemoryMiB          int         `json:"memory_mib"`
	TimeoutSeconds     int         `json:"timeout_seconds"`
	Files              []fileEntry `json:"files"`
}

type launchReceipt struct {
	Kind                           string      `json:"kind"`
	FreezeSHA256                   string      `json:"freeze_sha256"`
	SandboxID                      *string     `json:"sandbox_id"`
	StartedUnixNS                  int64       `json:"started_unix_ns"`
	FinishedUnixNS                 int64       `json:"finished_unix_ns"`
	Stage                          string      `json:"stage"`
	RemoteExitCode                 *int        `json:"remote_exit_code"`
	RemoteStdout                   string      `json:"remote_stdout"`
	RemoteStderr                   string      `json:"remote_stderr"`
	GuestReceiptPath               *string     `json:"guest_receipt_path"`
	GuestReceiptSHA                *string     `json:"guest_receipt_sha256"`
	Error                          *string     `json:"error"`
	TerminationError               *string     `json:"termination_error"`
	Files                          []fileEntry `json:"files"`
	HostExclusivePartitionVerified bool        `json:"host_exclusive_partition_verified"`
	ControlledSpeedup              *float64    `json:"controlled_speedup"`
}

func sha(raw []byte) string {
	h := sha256.Sum256(raw)
	return hex.EncodeToString(h[:])
}

func fileSHA(path string) (string, int64, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return "", 0, err
	}
	return sha(raw), int64(len(raw)), nil
}

func writeJSON(path string, value any) error {
	raw, err := json.MarshalIndent(value, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(path, append(raw, '\n'), 0600)
}

func preflight(here string) (plan, string, error) {
	var p plan
	freeze := filepath.Join(here, "PILOT_FREEZE.json")
	raw, err := os.ReadFile(freeze)
	if err != nil {
		return p, "", err
	}
	if err := json.Unmarshal(raw, &p); err != nil {
		return p, "", err
	}
	if p.Kind != "s3_modal_c1800_native_pilot_freeze_v1" || p.SDKVersion != "v0.11.0" ||
		p.GoVersion != runtime.Version() || p.SandboxCreationAttempts != 1 ||
		p.CPUPhysicalCores != 8 || p.MemoryMiB != 32768 ||
		p.LogicalExecutionCPUsMax != 16 || p.SandboxTimeoutSeconds != 1200 ||
		p.WorkloadDeadlineSeconds != 1080 || p.PerRunTimeoutSeconds != 240 ||
		p.MinimumRunSeconds != 30 || !p.BlockNetwork || !p.VMRuntime ||
		p.Region != nil || p.Image != "debian:bookworm-slim" ||
		p.AppName != "s3-c1800-n53-t205-native-pilot" || len(p.Files) != 4 {
		return p, "", errors.New("frozen resource or transfer contract mismatch")
	}
	checks := map[string]string{
		"go.mod": p.GoModSHA256, "go.sum": p.GoSumSHA256,
		"cmd/launch/main.go":                p.LaunchSourceSHA256,
		"cmd/analyze/main.go":               p.AnalyzerSourceSHA256,
		"BUILD_RECEIPT.json":                p.GuestBuildReceiptSHA256,
		"../modal/pilot-T205-20261010.json": p.PreviousReceiptSHA256,
	}
	for name, digest := range p.GuestSourceSHA256 {
		checks[name] = digest
	}
	for name, want := range checks {
		got, _, err := fileSHA(filepath.Join(here, name))
		if err != nil || got != want {
			return p, "", fmt.Errorf("source receipt mismatch %s: %s: %v", name, got, err)
		}
	}
	remoteNames := []string{"/opt/s3/baseline", "/opt/s3/expanded",
		"/opt/s3/public_bundle.json", "/opt/s3/guest"}
	var total int64
	for i, entry := range p.Files {
		if entry.Remote != remoteNames[i] {
			return p, "", fmt.Errorf("unexpected remote file %d: %s", i, entry.Remote)
		}
		local := entry.Local
		if !filepath.IsAbs(local) {
			local = filepath.Join(here, local)
		}
		got, size, err := fileSHA(local)
		if err != nil || got != entry.SHA256 || size != entry.Bytes {
			return p, "", fmt.Errorf("transfer mismatch %s: got %s/%d: %v", entry.Remote, got, size, err)
		}
		total += size
	}
	if total != p.TransferredBytes || p.Files[0].SHA256 != p.BinarySHA256["baseline"] ||
		p.Files[1].SHA256 != p.BinarySHA256["expanded"] ||
		p.Files[2].SHA256 != p.BundleSHA256 || p.Files[3].SHA256 != p.GuestSHA256 {
		return p, "", errors.New("transfer inventory hash or byte count mismatch")
	}
	bundle, err := os.ReadFile(filepath.Join(here, p.Files[2].Local))
	if err != nil || bytes.Contains(bundle, []byte("fixture_scalar")) {
		return p, "", errors.New("public bundle missing or contains fixture scalar")
	}
	var payload struct {
		Kind     string              `json:"kind"`
		RunOrder [][]json.RawMessage `json:"run_order"`
	}
	if json.Unmarshal(bundle, &payload) != nil || payload.Kind !=
		"s3_modal_c1800_public_t205_bundle_v1" || len(payload.RunOrder) != 4 {
		return p, "", errors.New("public bundle identity or ABBA order mismatch")
	}
	for i, want := range []string{"baseline", "expanded", "expanded", "baseline"} {
		var n int
		var arm string
		if len(payload.RunOrder[i]) != 2 || json.Unmarshal(payload.RunOrder[i][0], &n) != nil ||
			json.Unmarshal(payload.RunOrder[i][1], &arm) != nil || n != 53 || arm != want {
			return p, "", errors.New("public bundle ABBA order mismatch")
		}
	}
	return p, sha(raw), nil
}

func localPath(here, path string) string {
	if filepath.IsAbs(path) {
		return path
	}
	return filepath.Join(here, path)
}

func execAndRead(ctx context.Context, sb *modal.Sandbox, argv []string) (int, string, string, error) {
	process, err := sb.Exec(ctx, argv, nil)
	if err != nil {
		return -1, "", "", err
	}
	code, waitErr := process.Wait(ctx, nil)
	stdout, readOutErr := io.ReadAll(process.Stdout)
	stderr, readErrErr := io.ReadAll(process.Stderr)
	if readOutErr != nil || readErrErr != nil {
		return code, string(stdout), string(stderr),
			fmt.Errorf("stream read: stdout=%v stderr=%v", readOutErr, readErrErr)
	}
	return code, string(stdout), string(stderr), waitErr
}

func launch(here, output string, p plan, freezeSHA string) error {
	if output == "" {
		return errors.New("--launch requires --output")
	}
	var err error
	output, err = filepath.Abs(output)
	if err != nil {
		return err
	}
	attemptPath := strings.TrimSuffix(output, filepath.Ext(output)) + ".attempt.json"
	guestPath := strings.TrimSuffix(output, filepath.Ext(output)) + ".guest.json"
	statePath := strings.TrimSuffix(output, filepath.Ext(output)) + ".state.json"
	for _, path := range []string{output, attemptPath, guestPath, statePath} {
		if _, err := os.Stat(path); err == nil {
			return fmt.Errorf("a prior launch artifact already exists: %s", path)
		} else if !errors.Is(err, os.ErrNotExist) {
			return err
		}
	}
	if err := os.MkdirAll(filepath.Dir(output), 0700); err != nil {
		return err
	}
	ctx, cancel := context.WithTimeout(context.Background(), 1230*time.Second)
	defer cancel()
	mc, err := modal.NewClient()
	if err != nil {
		return err
	}
	app, err := mc.Apps.FromName(ctx, p.AppName, &modal.AppFromNameParams{CreateIfMissing: true})
	if err != nil {
		return err
	}
	image := mc.Images.FromRegistry(p.Image, nil)
	started := time.Now().UnixNano()
	a := attempt{Kind: "s3_modal_c1800_native_attempt_v1", FreezeSHA256: freezeSHA,
		StartedUnixNS: started, SandboxCreateCalls: 1,
		CPUPhysicalCores: p.CPUPhysicalCores, MemoryMiB: p.MemoryMiB,
		TimeoutSeconds: p.SandboxTimeoutSeconds, Files: p.Files}
	if err := writeJSON(attemptPath, a); err != nil {
		return err
	}
	receipt := launchReceipt{Kind: "s3_modal_c1800_native_launch_receipt_v1",
		FreezeSHA256: freezeSHA, StartedUnixNS: started, Stage: "sandbox_create",
		Files: p.Files}
	var sb *modal.Sandbox
	defer func() {
		if sb != nil {
			termCtx, termCancel := context.WithTimeout(context.Background(), 20*time.Second)
			defer termCancel()
			if _, termErr := sb.Terminate(termCtx, nil); termErr != nil {
				s := termErr.Error()
				receipt.TerminationError = &s
			}
		}
		receipt.FinishedUnixNS = time.Now().UnixNano()
		if saveErr := writeJSON(output, receipt); saveErr != nil {
			fmt.Fprintln(os.Stderr, "durable launch receipt:", saveErr)
		}
	}()
	sb, err = mc.Sandboxes.Create(ctx, app, image, &modal.SandboxCreateParams{
		CPU: float64(p.CPUPhysicalCores), CPULimit: float64(p.CPUPhysicalCores),
		MemoryMiB: p.MemoryMiB, MemoryLimitMiB: p.MemoryMiB,
		Runtime:      modal.SandboxRuntimeVM,
		Timeout:      time.Duration(p.SandboxTimeoutSeconds) * time.Second,
		BlockNetwork: true, Command: []string{"sleep", strconv.Itoa(p.SandboxTimeoutSeconds)},
	})
	if err != nil {
		s := "sandbox_create:" + err.Error()
		receipt.Error = &s
		return err
	}
	id := sb.SandboxID
	receipt.SandboxID = &id
	_ = writeJSON(statePath, receipt)
	fmt.Printf("{\"sandbox_id\":%q}\n", id)
	receipt.Stage = "upload"
	for _, entry := range p.Files {
		if err := sb.Filesystem.CopyFromLocal(ctx, localPath(here, entry.Local), entry.Remote, nil); err != nil {
			s := "upload:" + entry.Remote + ":" + err.Error()
			receipt.Error = &s
			return err
		}
	}
	receipt.Stage = "guest_exec"
	_ = writeJSON(statePath, receipt)
	code, _, stderr, err := execAndRead(ctx, sb, []string{"chmod", "0755", "/opt/s3/guest"})
	if err != nil || code != 0 {
		s := fmt.Sprintf("guest_chmod:code=%d:error=%v:stderr=%s", code, err, stderr)
		receipt.Error = &s
		return errors.New(s)
	}
	expected := map[string]any{
		"bundle_sha256": p.BundleSHA256, "binary_sha256": p.BinarySHA256,
		"guest_sha256":               p.GuestSHA256,
		"logical_execution_cpus_max": p.LogicalExecutionCPUsMax,
		"physical_core_limit":        p.CPUPhysicalCores,
		"workload_deadline_seconds":  p.WorkloadDeadlineSeconds,
		"per_run_timeout_seconds":    p.PerRunTimeoutSeconds,
		"minimum_run_seconds":        p.MinimumRunSeconds,
	}
	expectedJSON, _ := json.Marshal(expected)
	code, stdout, stderr, err := execAndRead(ctx, sb,
		[]string{"/opt/s3/guest", string(expectedJSON)})
	receipt.RemoteExitCode = &code
	receipt.RemoteStdout = clip(stdout, 4000)
	receipt.RemoteStderr = clip(stderr, 4000)
	if err != nil || code != 0 {
		s := fmt.Sprintf("guest_exec:code=%d:error=%v", code, err)
		receipt.Error = &s
	}
	receipt.Stage = "guest_receipt_copy"
	if copyErr := sb.Filesystem.CopyToLocal(ctx, "/tmp/s3-modal-native-receipt.json", guestPath, nil); copyErr != nil {
		s := "guest_receipt_copy:" + copyErr.Error()
		if receipt.Error != nil {
			s = *receipt.Error + ";" + s
		}
		receipt.Error = &s
		return copyErr
	}
	guestSHA, _, err := fileSHA(guestPath)
	if err != nil {
		return err
	}
	receipt.GuestReceiptPath = &guestPath
	receipt.GuestReceiptSHA = &guestSHA
	receipt.Stage = "complete"
	if receipt.Error != nil {
		return errors.New(*receipt.Error)
	}
	return nil
}

func clip(s string, maxBytes int) string {
	if len(s) > maxBytes {
		return s[len(s)-maxBytes:]
	}
	return s
}

func main() {
	preview := flag.Bool("preview", false, "verify the complete frozen transfer without creating a Sandbox")
	paidLaunch := flag.Bool("launch", false, "make the one approved Sandbox creation attempt")
	output := flag.String("output", "", "durable launch receipt path")
	flag.Parse()
	if *preview == *paidLaunch {
		fmt.Fprintln(os.Stderr, "choose exactly one of --preview or --launch")
		os.Exit(2)
	}
	here, err := os.Getwd()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	p, freezeSHA, err := preflight(here)
	if err != nil {
		fmt.Fprintln(os.Stderr, "preflight:", err)
		os.Exit(1)
	}
	if *preview {
		_ = json.NewEncoder(os.Stdout).Encode(map[string]any{
			"freeze_sha256": freezeSHA, "files": p.Files,
			"transferred_bytes":                      p.TransferredBytes,
			"cpu_physical_cores_request_and_limit":   p.CPUPhysicalCores,
			"memory_mib_request_and_limit":           p.MemoryMiB,
			"sandbox_timeout_seconds":                p.SandboxTimeoutSeconds,
			"published_base_compute_usd_for_timeout": p.PublishedBaseComputeUSD,
		})
		return
	}
	if err := launch(here, *output, p, freezeSHA); err != nil {
		fmt.Fprintln(os.Stderr, "launch:", err)
		os.Exit(1)
	}
}
