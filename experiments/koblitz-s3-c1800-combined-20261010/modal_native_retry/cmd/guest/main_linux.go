//go:build linux && amd64

package main

import (
	"bytes"
	"compress/gzip"
	"crypto/sha256"
	"encoding/base64"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"runtime"
	"slices"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/aburan28/crypto/s3-modal-native-retry/internal/topology"
)

const base = "/opt/s3"
const receiptPath = "/tmp/s3-modal-native-receipt.json"

var phases = []string{"target_query", "target_pdp_charged", "target_relation_check",
	"target_descent", "target_recovery_check"}

type expected struct {
	BundleSHA256            string            `json:"bundle_sha256"`
	BinarySHA256            map[string]string `json:"binary_sha256"`
	GuestSHA256             string            `json:"guest_sha256"`
	LogicalExecutionCPUsMax int               `json:"logical_execution_cpus_max"`
	PhysicalCoreLimit       int               `json:"physical_core_limit"`
	WorkloadDeadlineSeconds int               `json:"workload_deadline_seconds"`
	PerRunTimeoutSeconds    int               `json:"per_run_timeout_seconds"`
	MinimumRunSeconds       int               `json:"minimum_run_seconds"`
}

type arm struct {
	CandidateID               string `json:"candidate_id"`
	Columns                   int    `json:"columns"`
	FactorBaseDigest          string `json:"factor_base_digest"`
	FactorBasePoints          int    `json:"factor_base_points"`
	PrepareRank               int    `json:"prepare_rank"`
	RootTableSlots            int    `json:"root_table_slots"`
	TargetSearchPrefixColumns int    `json:"target_search_prefix_columns"`
	TargetSearchSchedule      string `json:"target_search_schedule"`
}

type curve struct {
	N           int            `json:"n"`
	Label       string         `json:"label"`
	PublicPoint []uint64       `json:"public_point"`
	WorkloadID  string         `json:"workload_id"`
	Arms        map[string]arm `json:"arms"`
}

type bundle struct {
	Kind            string              `json:"kind"`
	Curves          []curve             `json:"curves"`
	RelationSeed    int                 `json:"relation_seed"`
	RelationWorkers int                 `json:"relation_workers"`
	TargetWorkers   int                 `json:"target_workers"`
	RunOrder        [][]json.RawMessage `json:"run_order"`
}

type nativeReport struct {
	RecoveredScalar              uint64             `json:"recovered_scalar"`
	GroupVerified                bool               `json:"group_verified"`
	TargetCount                  int                `json:"target_count"`
	Target                       []uint64           `json:"target"`
	FactorBaseDigest             string             `json:"factor_base_digest"`
	FactorBasePoints             int                `json:"factor_base_points"`
	OrbitColumns                 int                `json:"orbit_columns"`
	PreparedRank                 int                `json:"prepared_rank"`
	RootTableSlots               int                `json:"root_table_slots"`
	TargetQueryWorkers           int                `json:"target_query_workers"`
	FallbackAttempts             int                `json:"fallback_attempts"`
	TargetFullRankDirectRecovery bool               `json:"target_full_rank_direct_recovery"`
	TargetSearchSchedule         string             `json:"target_search_schedule"`
	TargetSearchPrefixColumns    int                `json:"target_search_prefix_columns"`
	TimingMS                     map[string]float64 `json:"timing_ms"`
}

type host struct {
	Platform               string                     `json:"platform"`
	CPUModel               *string                    `json:"cpu_model"`
	CPUCountVisible        int                        `json:"cpu_count_visible"`
	AllowedCPUs            []int                      `json:"allowed_cpus"`
	SelectedCPUs           []int                      `json:"selected_cpus"`
	Selection              topology.Selection         `json:"selection"`
	Places                 map[int]topology.Place     `json:"places"`
	NUMANodes              map[int][]int              `json:"numa_nodes"`
	IRQRoutingTop64        map[string]irqRoute        `json:"irq_routing_top64"`
	CPUGovernor            map[int]*string            `json:"cpu_governor"`
	CPUFrequencyLimitsKHz  map[int]map[string]*string `json:"cpu_frequency_limits_khz"`
	ProcessMemoryAllowed   *string                    `json:"process_memory_allowed"`
	PerfEventParanoid      *string                    `json:"perf_event_paranoid"`
	HostExclusivePartition *string                    `json:"host_exclusive_partition_receipt"`
	IdleBefore             snapshot                   `json:"idle_before"`
	IdleAfter              snapshot                   `json:"idle_after"`
}

type runRow struct {
	N                       int               `json:"n"`
	Label                   string            `json:"label"`
	Arm                     string            `json:"arm"`
	Repetition              int               `json:"repetition"`
	CandidateID             string            `json:"candidate_id"`
	WorkloadID              string            `json:"workload_id"`
	RunID                   string            `json:"run_id"`
	Argv                    []string          `json:"argv"`
	Affinity                []int             `json:"affinity"`
	StartedUnixNS           int64             `json:"started_unix_ns"`
	FinishedUnixNS          int64             `json:"finished_unix_ns"`
	OuterWallNS             int64             `json:"outer_wall_ns"`
	ExitCode                int               `json:"exit_code"`
	Error                   *string           `json:"error"`
	Stderr                  string            `json:"stderr"`
	RawSHA256               *string           `json:"raw_sha256"`
	ReportGzipBase64        string            `json:"report_gzip_base64"`
	NativeVerified          bool              `json:"native_verified"`
	PhaseOK                 bool              `json:"phase_ok"`
	SampledPeakRSSBytes     uint64            `json:"sampled_peak_rss_bytes"`
	SampledNUMAPages        map[string]uint64 `json:"sampled_numa_pages"`
	SampledProcessSchedstat *string           `json:"sampled_process_schedstat"`
	MaxThreadsSeen          int               `json:"max_threads_seen"`
	ThreadAffinityMasksSeen []string          `json:"thread_affinity_masks_seen"`
	Perf                    perfPair          `json:"perf"`
	Before                  snapshot          `json:"before"`
	After                   snapshot          `json:"after"`
	Report                  json.RawMessage   `json:"report"`
}

type receipt struct {
	Kind         string              `json:"kind"`
	BundleSHA256 string              `json:"bundle_sha256"`
	BinarySHA256 map[string]string   `json:"binary_sha256"`
	GuestSHA256  string              `json:"guest_sha256"`
	RunOrder     [][]json.RawMessage `json:"run_order"`
	Rows         []runRow            `json:"rows"`
	StopReason   *string             `json:"stop_reason"`
	GuestError   *string             `json:"guest_error"`
	Host         *host               `json:"host"`
}

func digest(raw []byte) string {
	h := sha256.Sum256(raw)
	return hex.EncodeToString(h[:])
}

func fileDigest(path string) (string, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return "", err
	}
	return digest(raw), nil
}

func save(result *receipt) error {
	raw, err := json.Marshal(result)
	if err != nil {
		return err
	}
	tmp := receiptPath + ".tmp"
	if err := os.WriteFile(tmp, append(raw, '\n'), 0600); err != nil {
		return err
	}
	return os.Rename(tmp, receiptPath)
}

func stop(result *receipt, reason, detail string) {
	result.StopReason = &reason
	if detail != "" {
		result.GuestError = &detail
	}
	_ = save(result)
}

func cpuModel() *string {
	raw := read("/proc/cpuinfo", 1_000_000)
	if raw == nil {
		return nil
	}
	for _, line := range strings.Split(*raw, "\n") {
		if strings.HasPrefix(line, "model name") {
			_, value, ok := strings.Cut(line, ":")
			if ok {
				model := strings.TrimSpace(value)
				return &model
			}
		}
	}
	return nil
}

func memoryAllowed() *string {
	raw := read("/proc/self/status", 50_000)
	if raw == nil {
		return nil
	}
	for _, line := range strings.Split(*raw, "\n") {
		if strings.HasPrefix(line, "Mems_allowed_list:") {
			value := strings.TrimSpace(strings.TrimPrefix(line, "Mems_allowed_list:"))
			return &value
		}
	}
	return nil
}

func makeHost(allowed []int, places map[int]topology.Place, nodes map[int][]int,
	selection topology.Selection) *host {
	governor := map[int]*string{}
	limits := map[int]map[string]*string{}
	for _, cpu := range selection.CPUs {
		prefix := fmt.Sprintf("/sys/devices/system/cpu/cpu%d/cpufreq", cpu)
		governor[cpu] = trim(filepath.Join(prefix, "scaling_governor"))
		limits[cpu] = map[string]*string{
			"scaling_min_freq": trim(filepath.Join(prefix, "scaling_min_freq")),
			"scaling_max_freq": trim(filepath.Join(prefix, "scaling_max_freq")),
		}
	}
	return &host{Platform: runtime.GOOS + "/" + runtime.GOARCH,
		CPUModel: cpuModel(), CPUCountVisible: runtime.NumCPU(),
		AllowedCPUs: allowed, SelectedCPUs: selection.CPUs,
		Selection: selection, Places: places, NUMANodes: nodes,
		IRQRoutingTop64: irqRouting(selection.CPUs),
		CPUGovernor:     governor, CPUFrequencyLimitsKHz: limits,
		ProcessMemoryAllowed: memoryAllowed(),
		PerfEventParanoid:    trim("/proc/sys/kernel/perf_event_paranoid")}
}

func phaseOK(report nativeReport) bool {
	timing := report.TimingMS
	if timing == nil {
		return false
	}
	sum := 0.0
	for _, key := range phases {
		value, ok := timing[key]
		if !ok || value < 0 {
			return false
		}
		sum += value
	}
	total, ok := timing["target_online_after_reusable_setup"]
	return ok && total >= 0 && sum-total < .02 && total-sum < .02
}

func truncated(raw []byte, limit int) string {
	if len(raw) > limit {
		raw = raw[len(raw)-limit:]
	}
	return string(raw)
}

func runOne(b bundle, c curve, armName string, repetition int, cpus []int,
	seconds time.Duration) runRow {
	a := c.Arms[armName]
	row := runRow{N: c.N, Label: c.Label, Arm: armName, Repetition: repetition,
		CandidateID: a.CandidateID, WorkloadID: c.WorkloadID,
		RunID:    a.CandidateID + "W" + c.WorkloadID + "R" + strconv.Itoa(repetition),
		ExitCode: -1, SampledNUMAPages: map[string]uint64{}}
	pointRaw, _ := json.Marshal(c.PublicPoint)
	pointPath := fmt.Sprintf("/tmp/s3_n%d_T205.json", c.N)
	outputPath := fmt.Sprintf("/tmp/s3_n%d_%s_R%d.json", c.N, armName, repetition)
	if err := os.WriteFile(pointPath, append(pointRaw, '\n'), 0600); err != nil {
		s := "point_write:" + err.Error()
		row.Error = &s
		return row
	}
	binary := filepath.Join(base, armName)
	args := []string{strconv.Itoa(c.N), "0", strconv.Itoa(a.Columns),
		strconv.Itoa(b.RelationSeed), pointPath, outputPath,
		strconv.Itoa(b.RelationWorkers), strconv.Itoa(a.PrepareRank),
		strconv.Itoa(b.TargetWorkers)}
	row.Argv = append([]string{binary}, args...)
	self, err := os.Executable()
	if err != nil {
		s := "guest_executable:" + err.Error()
		row.Error = &s
		return row
	}
	childArgs := append([]string{"--child", strings.Trim(strings.ReplaceAll(fmt.Sprint(cpus), " ", ","), "[]"), binary}, args...)
	cmd := exec.Command(self, childArgs...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout, cmd.Stderr = &stdout, &stderr
	row.Before = makeSnapshot(cpus)
	row.StartedUnixNS = time.Now().UnixNano()
	if err := cmd.Start(); err != nil {
		s := "child_start:" + err.Error()
		row.Error = &s
		row.FinishedUnixNS = time.Now().UnixNano()
		row.OuterWallNS = row.FinishedUnixNS - row.StartedUnixNS
		row.After = makeSnapshot(cpus)
		return row
	}
	// The child pins itself before exec. Recheck until that handoff has occurred.
	for i := 0; i < 100; i++ {
		row.Affinity, _ = ownAffinity(cmd.Process.Pid)
		if reflect.DeepEqual(row.Affinity, slices.Sorted(slices.Values(slices.Clone(cpus)))) {
			break
		}
		time.Sleep(time.Millisecond)
	}
	cyclesFD, cyclesErr := openPerf(cmd.Process.Pid, 0)
	instructionsFD, instructionsErr := openPerf(cmd.Process.Pid, 1)
	waited := make(chan error, 1)
	go func() { waited <- cmd.Wait() }()
	ticker := time.NewTicker(250 * time.Millisecond)
	timer := time.NewTimer(seconds)
	defer ticker.Stop()
	defer timer.Stop()
	polls := 0
	threadMasksSeen := map[string]bool{}
	waitErr := error(nil)
	for active := true; active; {
		select {
		case waitErr = <-waited:
			active = false
		case <-timer.C:
			s := "run_timeout"
			row.Error = &s
			_ = cmd.Process.Kill()
			waitErr = <-waited
			active = false
		case <-ticker.C:
			polls++
			row.SampledPeakRSSBytes = max(row.SampledPeakRSSBytes, processRSS(cmd.Process.Pid))
			row.SampledProcessSchedstat = read(fmt.Sprintf("/proc/%d/schedstat", cmd.Process.Pid), 200)
			if polls%4 == 0 {
				count, masks := threadMasks(cmd.Process.Pid)
				row.MaxThreadsSeen = max(row.MaxThreadsSeen, count)
				for _, mask := range masks {
					threadMasksSeen[mask] = true
				}
			}
			if polls%8 == 0 {
				row.SampledNUMAPages = numaPages(cmd.Process.Pid)
			}
		}
	}
	for mask := range threadMasksSeen {
		row.ThreadAffinityMasksSeen = append(row.ThreadAffinityMasksSeen, mask)
	}
	slices.Sort(row.ThreadAffinityMasksSeen)
	row.FinishedUnixNS = time.Now().UnixNano()
	row.OuterWallNS = row.FinishedUnixNS - row.StartedUnixNS
	row.After = makeSnapshot(cpus)
	row.Perf = finishPerf(cyclesFD, instructionsFD, cyclesErr, instructionsErr)
	row.ExitCode = cmd.ProcessState.ExitCode()
	row.Stderr = truncated(stderr.Bytes(), 4000)
	if waitErr != nil && row.Error == nil {
		s := "child_exit:" + waitErr.Error()
		row.Error = &s
	}
	raw, err := os.ReadFile(outputPath)
	if err != nil {
		if row.Error == nil {
			s := "report_missing:" + err.Error()
			row.Error = &s
		}
		return row
	}
	hash := digest(raw)
	row.RawSHA256 = &hash
	row.Report = raw
	var compressed bytes.Buffer
	zipper := gzip.NewWriter(&compressed)
	if _, err := zipper.Write(raw); err != nil {
		s := "report_compression:" + err.Error()
		row.Error = &s
		_ = zipper.Close()
		return row
	}
	if err := zipper.Close(); err != nil {
		s := "report_compression:" + err.Error()
		row.Error = &s
		return row
	}
	row.ReportGzipBase64 = base64.StdEncoding.EncodeToString(compressed.Bytes())
	var report nativeReport
	if err := json.Unmarshal(raw, &report); err != nil {
		s := "report_invalid_json:" + err.Error()
		row.Error = &s
		return row
	}
	row.PhaseOK = phaseOK(report)
	row.NativeVerified = row.ExitCode == 0 && row.Error == nil &&
		bytes.Equal(stdout.Bytes(), raw) && row.PhaseOK && report.GroupVerified &&
		report.TargetCount == 1 && reflect.DeepEqual(report.Target, c.PublicPoint) &&
		report.FactorBaseDigest == a.FactorBaseDigest &&
		report.FactorBasePoints == a.FactorBasePoints && report.OrbitColumns == a.Columns &&
		report.PreparedRank == a.PrepareRank && report.RootTableSlots == a.RootTableSlots &&
		report.TargetQueryWorkers == 1 && report.FallbackAttempts == 0 &&
		report.TargetFullRankDirectRecovery &&
		report.TargetSearchSchedule == a.TargetSearchSchedule &&
		report.TargetSearchPrefixColumns == a.TargetSearchPrefixColumns
	return row
}

func execute(e expected) receipt {
	result := receipt{Kind: "s3_modal_c1800_native_guest_receipt_v1",
		BundleSHA256: e.BundleSHA256, BinarySHA256: e.BinarySHA256,
		GuestSHA256: e.GuestSHA256, Rows: []runRow{}}
	if e.PhysicalCoreLimit != 8 || e.LogicalExecutionCPUsMax != 16 ||
		e.WorkloadDeadlineSeconds != 1080 || e.PerRunTimeoutSeconds != 240 ||
		e.MinimumRunSeconds != 30 {
		stop(&result, "invalid_resource_freeze", "guest received unexpected resource limits")
		return result
	}
	self, err := os.Executable()
	if err != nil {
		stop(&result, "guest_identity_failure", err.Error())
		return result
	}
	if got, err := fileDigest(self); err != nil || got != e.GuestSHA256 {
		stop(&result, "guest_identity_mismatch", fmt.Sprintf("got=%s:error=%v", got, err))
		return result
	}
	for remote, want := range map[string]string{
		"public_bundle.json": e.BundleSHA256,
		"baseline":           e.BinarySHA256["baseline"],
		"expanded":           e.BinarySHA256["expanded"]} {
		got, err := fileDigest(filepath.Join(base, remote))
		if err != nil || got != want {
			stop(&result, "transfer_hash_mismatch", fmt.Sprintf("%s:got=%s:error=%v", remote, got, err))
			return result
		}
	}
	bundleRaw, err := os.ReadFile(filepath.Join(base, "public_bundle.json"))
	if err != nil {
		stop(&result, "bundle_read_failure", err.Error())
		return result
	}
	var b bundle
	if err := json.Unmarshal(bundleRaw, &b); err != nil ||
		b.Kind != "s3_modal_c1800_public_t205_bundle_v1" || len(b.Curves) != 1 {
		stop(&result, "invalid_bundle", fmt.Sprintf("%v", err))
		return result
	}
	result.RunOrder = b.RunOrder
	allowed, err := ownAffinity(0)
	if err != nil {
		stop(&result, "affinity_read_failure", err.Error())
		return result
	}
	places, nodes := readTopology(allowed)
	selection := topology.Select(allowed, places, e.PhysicalCoreLimit, e.LogicalExecutionCPUsMax)
	result.Host = makeHost(allowed, places, nodes, selection)
	result.Host.IdleBefore = makeSnapshot(selection.CPUs)
	time.Sleep(2 * time.Second)
	result.Host.IdleAfter = makeSnapshot(selection.CPUs)
	_ = save(&result)
	if len(selection.CPUs) == 0 {
		stop(&result, "no_eligible_cpus", "guest affinity is empty")
		return result
	}
	for _, armName := range []string{"baseline", "expanded"} {
		if err := os.Chmod(filepath.Join(base, armName), 0755); err != nil {
			stop(&result, "binary_chmod_failure", err.Error())
			return result
		}
	}
	deadline := time.Now().Add(time.Duration(e.WorkloadDeadlineSeconds) * time.Second)
	repetitions := map[string]int{}
	for _, item := range b.RunOrder {
		if len(item) != 2 {
			stop(&result, "invalid_run_order", "run-order item has wrong arity")
			return result
		}
		var n int
		var armName string
		if json.Unmarshal(item[0], &n) != nil || json.Unmarshal(item[1], &armName) != nil ||
			n != b.Curves[0].N || (armName != "baseline" && armName != "expanded") {
			stop(&result, "invalid_run_order", "run-order item does not match T205")
			return result
		}
		remaining := time.Until(deadline)
		if remaining < time.Duration(e.MinimumRunSeconds)*time.Second {
			stop(&result, "workload_deadline", "less than minimum run time remains")
			return result
		}
		repetitions[armName]++
		limit := min(time.Duration(e.PerRunTimeoutSeconds)*time.Second, remaining-5*time.Second)
		row := runOne(b, b.Curves[0], armName, repetitions[armName], selection.CPUs, limit)
		result.Rows = append(result.Rows, row)
		_ = save(&result)
		if !row.NativeVerified {
			stop(&result, "invalid_or_timeout", "native report verification failed")
			return result
		}
	}
	_ = save(&result)
	return result
}

func main() {
	if len(os.Args) >= 4 && os.Args[1] == "--child" {
		cpus := cpuList(os.Args[2])
		if err := setAffinity(0, cpus); err != nil {
			fmt.Fprintln(os.Stderr, "child affinity:", err)
			os.Exit(101)
		}
		if err := syscall.Exec(os.Args[3], os.Args[3:], os.Environ()); err != nil {
			fmt.Fprintln(os.Stderr, "child exec:", err)
			os.Exit(102)
		}
	}
	if len(os.Args) != 2 {
		fmt.Fprintln(os.Stderr, "usage: guest EXPECTED_JSON")
		os.Exit(2)
	}
	var e expected
	if err := json.Unmarshal([]byte(os.Args[1]), &e); err != nil {
		fmt.Fprintln(os.Stderr, "expected JSON:", err)
		os.Exit(2)
	}
	var result receipt
	defer func() {
		if recovered := recover(); recovered != nil {
			if raw, err := os.ReadFile(receiptPath); err == nil {
				_ = json.Unmarshal(raw, &result)
			}
			if result.Kind == "" {
				result.Kind = "s3_modal_c1800_native_guest_receipt_v1"
			}
			reason := "guest_panic"
			detail := fmt.Sprint(recovered)
			result.StopReason = &reason
			result.GuestError = &detail
			_ = save(&result)
		}
	}()
	result = execute(e)
	if err := save(&result); err != nil {
		fmt.Fprintln(os.Stderr, "receipt write:", err)
		os.Exit(1)
	}
	fmt.Printf("{\"receipt\":%q,\"rows\":%d}\n", receiptPath, len(result.Rows))
}
