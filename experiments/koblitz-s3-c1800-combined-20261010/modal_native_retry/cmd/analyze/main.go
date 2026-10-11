package main

import (
	"bytes"
	"compress/gzip"
	"crypto/sha256"
	"encoding/base64"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
)

var phaseKeys = []string{"target_query", "target_pdp_charged", "target_relation_check",
	"target_descent", "target_recovery_check"}

var qualityKeys = []string{
	"verified_scalar_and_fixture", "exclusive_phase_sum", "eight_physical_cores_selected",
	"whole_smt_siblings_selected", "observed_thread_affinity", "observed_worker_threads",
	"one_numa_node", "numa_pages_on_selected_node", "unlimited_cpu_quota",
	"zero_cpu_throttling", "zero_guest_steal", "zero_selected_cpu_interrupts",
	"cpu_psi_at_most_1_percent", "memory_io_psi_at_most_0_1_percent",
	"zero_oom_events", "fixed_cpu_frequency", "perf_counters_unmultiplexed",
}

type curveArm struct {
	CandidateID               string `json:"candidate_id"`
	FactorBaseDigest          string `json:"factor_base_digest"`
	FactorBasePoints          int    `json:"factor_base_points"`
	Columns                   int    `json:"columns"`
	PrepareRank               int    `json:"prepare_rank"`
	RootTableSlots            int    `json:"root_table_slots"`
	TargetSearchSchedule      string `json:"target_search_schedule"`
	TargetSearchPrefixColumns int    `json:"target_search_prefix_columns"`
}
type curve struct {
	N           int                 `json:"n"`
	Label       string              `json:"label"`
	PublicPoint []uint64            `json:"public_point"`
	WorkloadID  string              `json:"workload_id"`
	Arms        map[string]curveArm `json:"arms"`
}
type bundle struct {
	Curves   []curve             `json:"curves"`
	RunOrder [][]json.RawMessage `json:"run_order"`
}
type fixture struct {
	N                 int      `json:"n"`
	Label             string   `json:"label"`
	FixtureScalar     uint64   `json:"fixture_scalar"`
	PointFileRelative string   `json:"point_file_relative"`
	PointSHA256       string   `json:"point_sha256"`
	PublicPoint       []uint64 `json:"public_point"`
}
type parent struct {
	Targets []fixture `json:"targets"`
}
type launch struct {
	Kind             string  `json:"kind"`
	FreezeSHA256     string  `json:"freeze_sha256"`
	SandboxID        *string `json:"sandbox_id"`
	RemoteExitCode   *int    `json:"remote_exit_code"`
	GuestReceiptPath *string `json:"guest_receipt_path"`
	GuestReceiptSHA  *string `json:"guest_receipt_sha256"`
	Error            *string `json:"error"`
}
type freeze struct {
	Kind             string            `json:"kind"`
	BundleSHA256     string            `json:"bundle_sha256"`
	BinarySHA256     map[string]string `json:"binary_sha256"`
	GuestSHA256      string            `json:"guest_sha256"`
	CPUPhysicalCores int               `json:"cpu_physical_cores_request_and_limit"`
}
type cpuPlace struct {
	CPU      int   `json:"cpu"`
	Package  *int  `json:"package"`
	Core     *int  `json:"core"`
	Node     *int  `json:"node"`
	Siblings []int `json:"siblings"`
}
type selection struct {
	CPUs                 []int  `json:"cpus"`
	KnownPhysicalCores   int    `json:"known_physical_cores"`
	CoreIdentityComplete bool   `json:"core_identity_complete"`
	WholeSMT             bool   `json:"whole_smt"`
	SingleNUMANode       bool   `json:"single_numa_node"`
	Mode                 string `json:"mode"`
}
type cpuTicks struct {
	Steal uint64 `json:"steal"`
}
type snapshot struct {
	UnixNS     int64              `json:"unix_ns"`
	CPUTicks   map[int]cpuTicks   `json:"cpu_ticks"`
	Interrupts map[int]uint64     `json:"interrupts"`
	Cgroup     map[string]*string `json:"cgroup"`
	Pressure   map[string]*string `json:"pressure"`
}
type host struct {
	AllowedCPUs           []int                      `json:"allowed_cpus"`
	SelectedCPUs          []int                      `json:"selected_cpus"`
	Selection             selection                  `json:"selection"`
	Places                map[int]cpuPlace           `json:"places"`
	NUMANodes             map[int][]int              `json:"numa_nodes"`
	CPUGovernor           map[int]*string            `json:"cpu_governor"`
	CPUFrequencyLimitsKHz map[int]map[string]*string `json:"cpu_frequency_limits_khz"`
	IdleBefore            snapshot                   `json:"idle_before"`
	IdleAfter             snapshot                   `json:"idle_after"`
}
type perfEvent struct {
	RunningFrac float64 `json:"running_fraction"`
}
type perfPair struct {
	Cycles            *perfEvent `json:"cycles_event"`
	Instructions      *perfEvent `json:"instructions_event"`
	CyclesError       *string    `json:"cycles_error"`
	InstructionsError *string    `json:"instructions_error"`
	IPC               *float64   `json:"ipc"`
	IPCScaledEstimate *float64   `json:"ipc_scaled_estimate"`
	BothAtLeast99Pct  bool       `json:"both_events_at_least_99_percent_running"`
}
type row struct {
	N                       int               `json:"n"`
	Label                   string            `json:"label"`
	Arm                     string            `json:"arm"`
	Repetition              int               `json:"repetition"`
	CandidateID             string            `json:"candidate_id"`
	WorkloadID              string            `json:"workload_id"`
	RunID                   string            `json:"run_id"`
	Affinity                []int             `json:"affinity"`
	OuterWallNS             int64             `json:"outer_wall_ns"`
	ExitCode                int               `json:"exit_code"`
	Error                   *string           `json:"error"`
	RawSHA256               *string           `json:"raw_sha256"`
	ReportGzipBase64        string            `json:"report_gzip_base64"`
	NativeVerified          bool              `json:"native_verified"`
	PhaseOK                 bool              `json:"phase_ok"`
	SampledPeakRSSBytes     uint64            `json:"sampled_peak_rss_bytes"`
	SampledNUMAPages        map[string]uint64 `json:"sampled_numa_pages"`
	MaxThreadsSeen          int               `json:"max_threads_seen"`
	ThreadAffinityMasksSeen []string          `json:"thread_affinity_masks_seen"`
	Perf                    perfPair          `json:"perf"`
	Before                  snapshot          `json:"before"`
	After                   snapshot          `json:"after"`
	Report                  json.RawMessage   `json:"report"`
}
type guest struct {
	Kind         string              `json:"kind"`
	BundleSHA256 string              `json:"bundle_sha256"`
	BinarySHA256 map[string]string   `json:"binary_sha256"`
	GuestSHA256  string              `json:"guest_sha256"`
	RunOrder     [][]json.RawMessage `json:"run_order"`
	Rows         []row               `json:"rows"`
	StopReason   *string             `json:"stop_reason"`
	GuestError   *string             `json:"guest_error"`
	Host         *host               `json:"host"`
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
type quality struct {
	GuestScore0To100    float64          `json:"guest_score_0_to_100"`
	AvailableCheckCount int              `json:"available_check_count"`
	Checks              map[string]*bool `json:"checks"`
	FailedChecks        []string         `json:"failed_checks"`
	UnavailableChecks   []string         `json:"unavailable_checks"`
	Counters            map[string]any   `json:"counters"`
}
type rowSummary struct {
	N                        int             `json:"n"`
	Label                    string          `json:"label"`
	Arm                      string          `json:"arm"`
	RunID                    string          `json:"run_id"`
	VerifiedNativeAndFixture bool            `json:"verified_native_and_fixture"`
	ValidationChecks         map[string]bool `json:"validation_checks"`
	OnlineMS                 *float64        `json:"online_ms"`
	SetupMS                  *float64        `json:"setup_ms"`
	SampledPeakRSSBytes      uint64          `json:"sampled_peak_rss_bytes"`
	IPC                      *float64        `json:"ipc"`
	IPCScaledEstimate        *float64        `json:"ipc_scaled_estimate"`
	GuestQuality             quality         `json:"guest_quality"`
	Error                    *string         `json:"error"`
}
type curveSummary struct {
	BaselineOnlineMS             []float64          `json:"baseline_online_ms"`
	ExpandedOnlineMS             []float64          `json:"expanded_online_ms"`
	PairedExploratoryOnlineRatio float64            `json:"paired_exploratory_online_ratio"`
	AAMaxMinRatio                map[string]float64 `json:"aa_max_min_ratio"`
	AAWithin10Percent            bool               `json:"aa_within_10_percent"`
}
type analysis struct {
	Kind                           string                  `json:"kind"`
	LaunchReceiptSHA256            string                  `json:"launch_receipt_sha256"`
	Status                         string                  `json:"status"`
	Error                          *string                 `json:"error"`
	PlannedRows                    int                     `json:"planned_rows"`
	CompletedRows                  int                     `json:"completed_rows"`
	GuestStopReason                *string                 `json:"guest_stop_reason"`
	GuestError                     *string                 `json:"guest_error"`
	GuestAllowedCPUs               []int                   `json:"guest_allowed_cpus"`
	GuestSelectedCPUs              []int                   `json:"guest_selected_cpus"`
	GuestSelectionMode             *string                 `json:"guest_selection_mode"`
	IdleNoise                      map[string]any          `json:"idle_noise"`
	Rows                           []rowSummary            `json:"rows"`
	Curves                         map[string]curveSummary `json:"curves"`
	GuestScoreMin                  *float64                `json:"guest_score_min"`
	GuestScoreMean                 *float64                `json:"guest_score_mean"`
	HostExclusivePartitionVerified bool                    `json:"host_exclusive_partition_verified"`
	ControlledSpeedup              *float64                `json:"controlled_speedup"`
}

func sha(raw []byte) string {
	h := sha256.Sum256(raw)
	return hex.EncodeToString(h[:])
}
func load(path string, value any) ([]byte, error) {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	return raw, json.Unmarshal(raw, value)
}
func boolPtr(v bool) *bool { return &v }
func first(raw *string) string {
	if raw == nil {
		return ""
	}
	fields := strings.Fields(*raw)
	if len(fields) == 0 {
		return ""
	}
	return fields[0]
}
func named(raw *string, key string) *uint64 {
	if raw == nil {
		return nil
	}
	for _, line := range strings.Split(*raw, "\n") {
		fields := strings.Fields(line)
		if len(fields) >= 2 && fields[0] == key {
			if value, err := strconv.ParseUint(fields[1], 10, 64); err == nil {
				return &value
			}
		}
	}
	return nil
}
func delta(before, after *string, key string) *uint64 {
	a, b := named(before, key), named(after, key)
	if a == nil || b == nil || *b < *a {
		return nil
	}
	v := *b - *a
	return &v
}
func pressureTotal(raw *string) *uint64 {
	if raw == nil {
		return nil
	}
	for _, line := range strings.Split(*raw, "\n") {
		if !strings.HasPrefix(line, "some ") {
			continue
		}
		for _, field := range strings.Fields(line) {
			if value, ok := strings.CutPrefix(field, "total="); ok {
				if n, err := strconv.ParseUint(value, 10, 64); err == nil {
					return &n
				}
			}
		}
	}
	return nil
}
func psiFraction(before, after *string, outerNS int64) *float64 {
	a, b := pressureTotal(before), pressureTotal(after)
	if a == nil || b == nil || *b < *a || outerNS <= 0 {
		return nil
	}
	v := float64(*b-*a) / (float64(outerNS) / 1000)
	return &v
}
func cpuList(raw string) []int {
	seen := map[int]bool{}
	for _, term := range strings.Split(strings.TrimSpace(raw), ",") {
		if term == "" {
			continue
		}
		lo, hi := term, term
		if a, b, ok := strings.Cut(term, "-"); ok {
			lo, hi = a, b
		}
		a, e1 := strconv.Atoi(strings.TrimSpace(lo))
		b, e2 := strconv.Atoi(strings.TrimSpace(hi))
		if e1 != nil || e2 != nil || a < 0 || b < a || b > 4096 {
			continue
		}
		for cpu := a; cpu <= b; cpu++ {
			seen[cpu] = true
		}
	}
	out := make([]int, 0, len(seen))
	for cpu := range seen {
		out = append(out, cpu)
	}
	sort.Ints(out)
	return out
}
func selectedTickDelta(before, after snapshot, cpus []int) *uint64 {
	var total uint64
	for _, cpu := range cpus {
		a, okA := before.CPUTicks[cpu]
		b, okB := after.CPUTicks[cpu]
		if !okA || !okB || b.Steal < a.Steal {
			return nil
		}
		total += b.Steal - a.Steal
	}
	return &total
}
func selectedIRQDelta(before, after snapshot, cpus []int) *uint64 {
	if before.Interrupts == nil || after.Interrupts == nil {
		return nil
	}
	var total uint64
	for _, cpu := range cpus {
		a, okA := before.Interrupts[cpu]
		b, okB := after.Interrupts[cpu]
		if !okA || !okB || b < a {
			return nil
		}
		total += b - a
	}
	return &total
}
func phaseOK(report nativeReport) bool {
	if report.TimingMS == nil {
		return false
	}
	sum := 0.0
	for _, key := range phaseKeys {
		v, ok := report.TimingMS[key]
		if !ok || v < 0 {
			return false
		}
		sum += v
	}
	v, ok := report.TimingMS["target_online_after_reusable_setup"]
	return ok && v > 0 && math.Abs(sum-v) < .02
}

func rawReportOK(r row) bool {
	if r.RawSHA256 == nil || r.ReportGzipBase64 == "" {
		return false
	}
	packed, err := base64.StdEncoding.DecodeString(r.ReportGzipBase64)
	if err != nil {
		return false
	}
	reader, err := gzip.NewReader(bytes.NewReader(packed))
	if err != nil {
		return false
	}
	defer reader.Close()
	raw, err := io.ReadAll(io.LimitReader(reader, 16*1024*1024+1))
	if err != nil || len(raw) > 16*1024*1024 || sha(raw) != *r.RawSHA256 {
		return false
	}
	var original, embedded bytes.Buffer
	return json.Compact(&original, raw) == nil && json.Compact(&embedded, r.Report) == nil &&
		bytes.Equal(original.Bytes(), embedded.Bytes())
}
func timing(report nativeReport, key string) *float64 {
	v, ok := report.TimingMS[key]
	if !ok {
		return nil
	}
	return &v
}
func score(row row, h host, verified bool, physicalRequested int) quality {
	selected := slices.Clone(h.SelectedCPUs)
	sort.Ints(selected)
	checks := map[string]*bool{}
	checks["verified_scalar_and_fixture"] = boolPtr(verified)
	checks["exclusive_phase_sum"] = boolPtr(row.PhaseOK)
	checks["eight_physical_cores_selected"] = boolPtr(
		h.Selection.CoreIdentityComplete && h.Selection.KnownPhysicalCores == physicalRequested)
	checks["whole_smt_siblings_selected"] = boolPtr(h.Selection.WholeSMT)
	if len(row.ThreadAffinityMasksSeen) > 0 {
		match := slices.Equal(row.Affinity, selected)
		for _, mask := range row.ThreadAffinityMasksSeen {
			match = match && slices.Equal(cpuList(mask), selected)
		}
		checks["observed_thread_affinity"] = boolPtr(match)
	}
	checks["observed_worker_threads"] = boolPtr(row.MaxThreadsSeen >= 14)
	if len(h.NUMANodes) > 0 {
		checks["one_numa_node"] = boolPtr(h.Selection.SingleNUMANode)
	}
	if len(row.SampledNUMAPages) > 0 && h.Selection.SingleNUMANode && len(selected) > 0 {
		p, ok := h.Places[selected[0]]
		if ok && p.Node != nil {
			want := "N" + strconv.Itoa(*p.Node)
			match := len(row.SampledNUMAPages) == 1
			for node := range row.SampledNUMAPages {
				match = match && node == want
			}
			checks["numa_pages_on_selected_node"] = boolPtr(match)
		}
	}
	if raw := row.Before.Cgroup["cpu.max"]; raw != nil {
		checks["unlimited_cpu_quota"] = boolPtr(first(raw) == "max")
	}
	throttled := delta(row.Before.Cgroup["cpu.stat"], row.After.Cgroup["cpu.stat"], "nr_throttled")
	throttleUS := delta(row.Before.Cgroup["cpu.stat"], row.After.Cgroup["cpu.stat"], "throttled_usec")
	if throttled != nil && throttleUS != nil {
		checks["zero_cpu_throttling"] = boolPtr(*throttled == 0 && *throttleUS == 0)
	}
	steal := selectedTickDelta(row.Before, row.After, selected)
	if steal != nil {
		checks["zero_guest_steal"] = boolPtr(*steal == 0)
	}
	irq := selectedIRQDelta(row.Before, row.After, selected)
	if irq != nil {
		checks["zero_selected_cpu_interrupts"] = boolPtr(*irq == 0)
	}
	cpuPSI := psiFraction(row.Before.Cgroup["cpu.pressure"], row.After.Cgroup["cpu.pressure"], row.OuterWallNS)
	memPSI := psiFraction(row.Before.Cgroup["memory.pressure"], row.After.Cgroup["memory.pressure"], row.OuterWallNS)
	ioPSI := psiFraction(row.Before.Cgroup["io.pressure"], row.After.Cgroup["io.pressure"], row.OuterWallNS)
	if cpuPSI != nil {
		checks["cpu_psi_at_most_1_percent"] = boolPtr(*cpuPSI >= 0 && *cpuPSI <= .01)
	}
	if memPSI != nil && ioPSI != nil {
		checks["memory_io_psi_at_most_0_1_percent"] = boolPtr(
			*memPSI >= 0 && *memPSI <= .001 && *ioPSI >= 0 && *ioPSI <= .001)
	}
	oom := delta(row.Before.Cgroup["memory.events"], row.After.Cgroup["memory.events"], "oom")
	oomKill := delta(row.Before.Cgroup["memory.events"], row.After.Cgroup["memory.events"], "oom_kill")
	if oom != nil && oomKill != nil {
		checks["zero_oom_events"] = boolPtr(*oom == 0 && *oomKill == 0)
	}
	freqKnown := len(selected) > 0
	fixed := true
	for _, cpu := range selected {
		governor := h.CPUGovernor[cpu]
		limits := h.CPUFrequencyLimitsKHz[cpu]
		low, high := limits["scaling_min_freq"], limits["scaling_max_freq"]
		if governor == nil || low == nil || high == nil {
			freqKnown = false
			break
		}
		fixed = fixed && strings.TrimSpace(*governor) == "performance" &&
			strings.TrimSpace(*low) == strings.TrimSpace(*high)
	}
	if freqKnown {
		checks["fixed_cpu_frequency"] = boolPtr(fixed)
	}
	if row.Perf.CyclesError == nil && row.Perf.InstructionsError == nil &&
		row.Perf.Cycles != nil && row.Perf.Instructions != nil {
		checks["perf_counters_unmultiplexed"] = boolPtr(row.Perf.BothAtLeast99Pct &&
			row.Perf.Cycles.RunningFrac >= .99 && row.Perf.Instructions.RunningFrac >= .99)
	}
	q := quality{Checks: checks, Counters: map[string]any{
		"cpu_throttled_events": throttled, "cpu_throttled_us": throttleUS,
		"selected_cpu_steal_ticks": steal, "selected_cpu_interrupts": irq,
		"cgroup_cpu_psi_fraction": cpuPSI, "cgroup_memory_psi_fraction": memPSI,
		"cgroup_io_psi_fraction": ioPSI, "oom_events": oom, "oom_kills": oomKill}}
	trueCount := 0
	for _, key := range qualityKeys {
		v := checks[key]
		if v == nil {
			q.UnavailableChecks = append(q.UnavailableChecks, key)
			continue
		}
		q.AvailableCheckCount++
		if *v {
			trueCount++
		} else {
			q.FailedChecks = append(q.FailedChecks, key)
		}
	}
	q.GuestScore0To100 = math.Round(10000*float64(trueCount)/float64(len(qualityKeys))) / 100
	return q
}

func analyze(here, launchPath string) (analysis, error) {
	result := analysis{Kind: "s3_modal_c1800_native_pilot_analysis_v1", Status: "producer_failure",
		PlannedRows: 4, Rows: []rowSummary{}, Curves: map[string]curveSummary{}}
	var p freeze
	freezeRaw, err := load(filepath.Join(here, "PILOT_FREEZE.json"), &p)
	if err != nil {
		return result, err
	}
	var l launch
	launchRaw, err := load(launchPath, &l)
	if err != nil {
		return result, err
	}
	result.LaunchReceiptSHA256 = sha(launchRaw)
	if l.Kind != "s3_modal_c1800_native_launch_receipt_v1" || l.FreezeSHA256 != sha(freezeRaw) {
		return result, errors.New("launch receipt does not match native freeze")
	}
	result.Error = l.Error
	if l.GuestReceiptPath == nil || l.GuestReceiptSHA == nil {
		return result, nil
	}
	var g guest
	guestRaw, err := load(*l.GuestReceiptPath, &g)
	if err != nil {
		return result, err
	}
	if sha(guestRaw) != *l.GuestReceiptSHA || g.Kind != "s3_modal_c1800_native_guest_receipt_v1" ||
		g.BundleSHA256 != p.BundleSHA256 || g.GuestSHA256 != p.GuestSHA256 ||
		g.BinarySHA256["baseline"] != p.BinarySHA256["baseline"] ||
		g.BinarySHA256["expanded"] != p.BinarySHA256["expanded"] {
		return result, errors.New("guest receipt identity mismatch")
	}
	var b bundle
	bundleRaw, err := load(filepath.Join(here, "../modal/public_bundle.json"), &b)
	if err != nil {
		return result, err
	}
	if sha(bundleRaw) != p.BundleSHA256 || len(b.Curves) != 1 || len(b.RunOrder) != 4 || len(g.RunOrder) != 4 {
		return result, errors.New("public bundle identity or run order mismatch")
	}
	for i := range b.RunOrder {
		if string(b.RunOrder[i][0]) != string(g.RunOrder[i][0]) ||
			string(b.RunOrder[i][1]) != string(g.RunOrder[i][1]) {
			return result, errors.New("guest run order mismatch")
		}
	}
	var parentFreeze parent
	_, err = load(filepath.Join(here, "../CONFIRM_FREEZE.json"), &parentFreeze)
	if err != nil {
		return result, err
	}
	var f *fixture
	for i := range parentFreeze.Targets {
		if parentFreeze.Targets[i].N == 53 && parentFreeze.Targets[i].Label == "T205" {
			f = &parentFreeze.Targets[i]
			break
		}
	}
	if f == nil {
		return result, errors.New("T205 fixture absent")
	}
	pointPath := filepath.Join(here, "../../..", f.PointFileRelative)
	pointRaw, err := os.ReadFile(pointPath)
	if err != nil {
		return result, err
	}
	var point []uint64
	pointOK := sha(pointRaw) == f.PointSHA256 && json.Unmarshal(pointRaw, &point) == nil &&
		slices.Equal(point, b.Curves[0].PublicPoint) && slices.Equal(point, f.PublicPoint)
	result.CompletedRows = len(g.Rows)
	result.GuestStopReason = g.StopReason
	result.GuestError = g.GuestError
	if g.Host != nil {
		result.GuestAllowedCPUs = g.Host.AllowedCPUs
		result.GuestSelectedCPUs = g.Host.SelectedCPUs
		mode := g.Host.Selection.Mode
		result.GuestSelectionMode = &mode
		result.IdleNoise = map[string]any{
			"elapsed_ns": g.Host.IdleAfter.UnixNS - g.Host.IdleBefore.UnixNS,
			"cpu_psi_fraction": psiFraction(g.Host.IdleBefore.Cgroup["cpu.pressure"],
				g.Host.IdleAfter.Cgroup["cpu.pressure"], g.Host.IdleAfter.UnixNS-g.Host.IdleBefore.UnixNS),
			"selected_cpu_steal_ticks": selectedTickDelta(g.Host.IdleBefore, g.Host.IdleAfter, g.Host.SelectedCPUs),
			"selected_cpu_interrupts":  selectedIRQDelta(g.Host.IdleBefore, g.Host.IdleAfter, g.Host.SelectedCPUs)}
	}
	verifiedAll := len(g.Rows) == 4 && l.Error == nil && l.SandboxID != nil &&
		l.RemoteExitCode != nil && *l.RemoteExitCode == 0 && g.Host != nil
	repetitions := map[string]int{}
	for i, r := range g.Rows {
		if i >= 4 {
			verifiedAll = false
			break
		}
		var wantN int
		var wantArm string
		if len(b.RunOrder[i]) != 2 || json.Unmarshal(b.RunOrder[i][0], &wantN) != nil ||
			json.Unmarshal(b.RunOrder[i][1], &wantArm) != nil {
			return result, errors.New("invalid frozen run order")
		}
		c := b.Curves[0]
		a := c.Arms[wantArm]
		repetitions[wantArm]++
		var report nativeReport
		reportOK := json.Unmarshal(r.Report, &report) == nil
		sumsOK := reportOK && phaseOK(report)
		rawOK := rawReportOK(r)
		validation := map[string]bool{
			"public_point_fixture":      pointOK,
			"guest_native_verification": r.NativeVerified,
			"exclusive_phase_sum":       r.PhaseOK && sumsOK,
			"raw_report_hash":           rawOK,
			"row_process_success":       r.Error == nil && r.ExitCode == 0 && r.OuterWallNS > 0,
			"run_identity": r.N == wantN && r.Arm == wantArm && r.Label == "T205" &&
				r.Repetition == repetitions[wantArm] &&
				r.CandidateID == a.CandidateID && r.WorkloadID == c.WorkloadID &&
				r.RunID == a.CandidateID+"W"+c.WorkloadID+"R"+strconv.Itoa(r.Repetition),
			"scalar_fixture": reportOK && report.RecoveredScalar == f.FixtureScalar,
			"point_replay": reportOK && report.GroupVerified && report.TargetCount == 1 &&
				slices.Equal(report.Target, c.PublicPoint),
			"factor_base_and_rank": reportOK && report.FactorBaseDigest == a.FactorBaseDigest &&
				report.FactorBasePoints == a.FactorBasePoints && report.OrbitColumns == a.Columns &&
				report.PreparedRank == a.PrepareRank && report.RootTableSlots == a.RootTableSlots,
			"target_search_contract": reportOK && report.TargetQueryWorkers == 1 &&
				report.FallbackAttempts == 0 && report.TargetFullRankDirectRecovery &&
				report.TargetSearchSchedule == a.TargetSearchSchedule &&
				report.TargetSearchPrefixColumns == a.TargetSearchPrefixColumns,
		}
		verified := true
		for _, ok := range validation {
			verified = verified && ok
		}
		if !verified {
			verifiedAll = false
		}
		var q quality
		if g.Host != nil {
			q = score(r, *g.Host, verified, p.CPUPhysicalCores)
		}
		result.Rows = append(result.Rows, rowSummary{N: r.N, Label: r.Label, Arm: r.Arm, RunID: r.RunID,
			VerifiedNativeAndFixture: verified, ValidationChecks: validation,
			OnlineMS: timing(report, "target_online_after_reusable_setup"),
			SetupMS:  timing(report, "reusable_setup_total"), SampledPeakRSSBytes: r.SampledPeakRSSBytes,
			IPC: r.Perf.IPC, IPCScaledEstimate: r.Perf.IPCScaledEstimate, GuestQuality: q, Error: r.Error})
	}
	if verifiedAll {
		arms := map[string][]float64{"baseline": {}, "expanded": {}}
		for _, r := range result.Rows {
			arms[r.Arm] = append(arms[r.Arm], *r.OnlineMS)
		}
		if len(arms["baseline"]) != 2 || len(arms["expanded"]) != 2 {
			return result, errors.New("paired arm count mismatch")
		}
		for _, times := range arms {
			for _, v := range times {
				if v <= 0 {
					return result, errors.New("invalid online timing")
				}
			}
		}
		aa := map[string]float64{}
		within := true
		for arm, times := range arms {
			aa[arm] = math.Max(times[0], times[1]) / math.Min(times[0], times[1])
			within = within && aa[arm] <= 1.10
		}
		result.Curves["53"] = curveSummary{BaselineOnlineMS: arms["baseline"], ExpandedOnlineMS: arms["expanded"],
			PairedExploratoryOnlineRatio: math.Sqrt(arms["baseline"][0]*arms["baseline"][1]) /
				math.Sqrt(arms["expanded"][0]*arms["expanded"][1]), AAMaxMinRatio: aa, AAWithin10Percent: within}
		result.Status = "verified_complete"
	} else if len(g.Rows) > 0 {
		result.Status = "incomplete_or_invalid"
	} else {
		result.Status = "no_rows"
	}
	if len(result.Rows) > 0 {
		minScore := result.Rows[0].GuestQuality.GuestScore0To100
		sum := 0.0
		for _, r := range result.Rows {
			minScore = math.Min(minScore, r.GuestQuality.GuestScore0To100)
			sum += r.GuestQuality.GuestScore0To100
		}
		mean := sum / float64(len(result.Rows))
		result.GuestScoreMin = &minScore
		result.GuestScoreMean = &mean
	}
	return result, nil
}

func main() {
	output := flag.String("output", "", "analysis artifact path")
	flag.Parse()
	if flag.NArg() != 1 || *output == "" {
		fmt.Fprintln(os.Stderr, "usage: analyze --output RESULT.json LAUNCH.json")
		os.Exit(2)
	}
	if _, err := os.Stat(*output); err == nil {
		fmt.Fprintln(os.Stderr, "analysis already exists")
		os.Exit(1)
	}
	here, err := os.Getwd()
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	result, err := analyze(here, flag.Arg(0))
	if err != nil {
		fmt.Fprintln(os.Stderr, "analyze:", err)
		os.Exit(1)
	}
	raw, err := json.MarshalIndent(result, "", "  ")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	if err := os.WriteFile(*output, append(raw, '\n'), 0600); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	_ = json.NewEncoder(os.Stdout).Encode(map[string]any{"status": result.Status, "completed_rows": result.CompletedRows,
		"curves": result.Curves, "controlled_speedup": nil})
}
