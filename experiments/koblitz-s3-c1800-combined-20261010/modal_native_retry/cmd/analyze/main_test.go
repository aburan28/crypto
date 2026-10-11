package main

import (
	"bytes"
	"compress/gzip"
	"encoding/base64"
	"encoding/json"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"testing"
)

func text(s string) *string { return &s }
func integer(n int) *int    { return &n }

func packedReport(t *testing.T, raw []byte) string {
	t.Helper()
	var compressed bytes.Buffer
	writer := gzip.NewWriter(&compressed)
	if _, err := writer.Write(raw); err != nil {
		t.Fatal(err)
	}
	if err := writer.Close(); err != nil {
		t.Fatal(err)
	}
	return base64.StdEncoding.EncodeToString(compressed.Bytes())
}

func archivedReport(t *testing.T, here, armName string, repetition int) json.RawMessage {
	t.Helper()
	localName := armName
	if armName == "expanded" {
		localName = "combined"
	}
	path := filepath.Join(here, "../local/paired_confirm",
		"T205_"+localName+"_R"+strconv.Itoa(repetition), "report.json.gz")
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	reader, err := gzip.NewReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()
	decoded, err := io.ReadAll(reader)
	if err != nil {
		t.Fatal(err)
	}
	return decoded
}

func syntheticReceipt(t *testing.T, here string) (guest, launch) {
	t.Helper()
	var p freeze
	freezeRaw, err := load(filepath.Join(here, "PILOT_FREEZE.json"), &p)
	if err != nil {
		t.Fatal(err)
	}
	var b bundle
	if _, err := load(filepath.Join(here, "../modal/public_bundle.json"), &b); err != nil {
		t.Fatal(err)
	}
	selected := []int{0, 1, 2, 3, 4, 5, 6, 7}
	places := map[int]cpuPlace{}
	governor := map[int]*string{}
	limits := map[int]map[string]*string{}
	ticks := map[int]cpuTicks{}
	interrupts := map[int]uint64{}
	for _, cpu := range selected {
		places[cpu] = cpuPlace{CPU: cpu, Package: integer(0), Core: integer(cpu),
			Node: integer(0), Siblings: []int{cpu, cpu + 8}}
		governor[cpu] = text("performance")
		limits[cpu] = map[string]*string{"scaling_min_freq": text("2500000"),
			"scaling_max_freq": text("2500000")}
		ticks[cpu] = cpuTicks{Steal: 0}
		interrupts[cpu] = 0
	}
	cgroup := map[string]*string{
		"cpu.max": text("max 100000"), "cpu.stat": text("nr_throttled 0\nthrottled_usec 0"),
		"cpu.pressure":    text("some avg10=0.00 total=0"),
		"memory.pressure": text("some avg10=0.00 total=0"),
		"io.pressure":     text("some avg10=0.00 total=0"),
		"memory.events":   text("oom 0\noom_kill 0"),
	}
	before := snapshot{UnixNS: 1_000_000_000, CPUTicks: ticks, Interrupts: interrupts,
		Cgroup: cgroup}
	after := before
	after.UnixNS += 30_000_000_000
	h := &host{AllowedCPUs: selected, SelectedCPUs: selected,
		Selection: selection{CPUs: selected, KnownPhysicalCores: 8,
			CoreIdentityComplete: true, WholeSMT: false,
			SingleNUMANode: true, Mode: "observed_affinity_cores"},
		Places: places, NUMANodes: map[int][]int{0: selected},
		CPUGovernor: governor, CPUFrequencyLimitsKHz: limits,
		IdleBefore: before, IdleAfter: after}
	g := guest{Kind: "s3_modal_c1800_native_guest_receipt_v1",
		BundleSHA256: p.BundleSHA256, BinarySHA256: p.BinarySHA256,
		GuestSHA256: p.GuestSHA256, RunOrder: b.RunOrder, Host: h}
	repetitions := map[string]int{}
	for _, item := range b.RunOrder {
		var armName string
		if err := json.Unmarshal(item[1], &armName); err != nil {
			t.Fatal(err)
		}
		repetitions[armName]++
		repetition := repetitions[armName]
		arm := b.Curves[0].Arms[armName]
		raw := archivedReport(t, here, armName, repetition)
		digest := sha(raw)
		g.Rows = append(g.Rows, row{N: 53, Label: "T205", Arm: armName,
			Repetition: repetition, CandidateID: arm.CandidateID,
			WorkloadID: b.Curves[0].WorkloadID,
			RunID: arm.CandidateID + "W" + b.Curves[0].WorkloadID +
				"R" + strconv.Itoa(repetition), Affinity: selected,
			OuterWallNS: 30_000_000_000, NativeVerified: true, PhaseOK: true,
			RawSHA256: &digest, Report: raw, ReportGzipBase64: packedReport(t, raw),
			SampledNUMAPages: map[string]uint64{"N0": 1000},
			MaxThreadsSeen:   14, ThreadAffinityMasksSeen: []string{"0-7"},
			Before: before, After: after, Perf: perfPair{
				Cycles: &perfEvent{RunningFrac: 1}, Instructions: &perfEvent{RunningFrac: 1},
				BothAtLeast99Pct: true}})
	}
	id := "sb-synthetic"
	exit := 0
	l := launch{Kind: "s3_modal_c1800_native_launch_receipt_v1",
		FreezeSHA256: sha(freezeRaw), SandboxID: &id, RemoteExitCode: &exit}
	return g, l
}

func callAnalyzer(t *testing.T, here string, g guest, l launch) analysis {
	t.Helper()
	dir := t.TempDir()
	guestRaw, err := json.Marshal(g)
	if err != nil {
		t.Fatal(err)
	}
	guestPath := filepath.Join(dir, "guest.json")
	if err := os.WriteFile(guestPath, guestRaw, 0600); err != nil {
		t.Fatal(err)
	}
	digest := sha(guestRaw)
	l.GuestReceiptPath = &guestPath
	l.GuestReceiptSHA = &digest
	launchRaw, err := json.Marshal(l)
	if err != nil {
		t.Fatal(err)
	}
	launchPath := filepath.Join(dir, "launch.json")
	if err := os.WriteFile(launchPath, launchRaw, 0600); err != nil {
		t.Fatal(err)
	}
	result, err := analyze(here, launchPath)
	if err != nil {
		t.Fatal(err)
	}
	return result
}

func TestArchivedT205ReportsAndPartialSMTTopology(t *testing.T) {
	here := filepath.Clean("../..")
	g, l := syntheticReceipt(t, here)
	result := callAnalyzer(t, here, g, l)
	if result.Status != "verified_complete" || result.CompletedRows != 4 ||
		len(result.Curves["53"].BaselineOnlineMS) != 2 || result.ControlledSpeedup != nil {
		t.Fatalf("archive analysis failed: status=%s first_row_checks=%v",
			result.Status, result.Rows[0].ValidationChecks)
	}
	for _, r := range result.Rows {
		if !r.VerifiedNativeAndFixture || !*r.GuestQuality.Checks["eight_physical_cores_selected"] ||
			*r.GuestQuality.Checks["whole_smt_siblings_selected"] {
			t.Fatalf("partial SMT quality mismatch: %+v", r)
		}
	}
}

func TestScalarMismatchInvalidatesPairedResult(t *testing.T) {
	here := filepath.Clean("../..")
	g, l := syntheticReceipt(t, here)
	var report map[string]any
	if err := json.Unmarshal(g.Rows[0].Report, &report); err != nil {
		t.Fatal(err)
	}
	report["recovered_scalar"] = float64(1)
	raw, _ := json.Marshal(report)
	g.Rows[0].Report = raw
	digest := sha(raw)
	g.Rows[0].RawSHA256 = &digest
	g.Rows[0].ReportGzipBase64 = packedReport(t, raw)
	result := callAnalyzer(t, here, g, l)
	if result.Status != "incomplete_or_invalid" || result.Rows[0].VerifiedNativeAndFixture ||
		len(result.Curves) != 0 {
		t.Fatalf("scalar mismatch accepted: %+v", result)
	}
}

func TestReportHashMismatchInvalidatesPairedResult(t *testing.T) {
	here := filepath.Clean("../..")
	g, l := syntheticReceipt(t, here)
	g.Rows[1].RawSHA256 = text("wrong")
	result := callAnalyzer(t, here, g, l)
	if result.Status != "incomplete_or_invalid" || result.Rows[1].VerifiedNativeAndFixture {
		t.Fatalf("raw report mismatch accepted: %+v", result)
	}
}
