//go:build linux && amd64

package main

import (
	"encoding/binary"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
	"syscall"
	"time"
	"unsafe"

	"github.com/aburan28/crypto/s3-modal-native-retry/internal/topology"
)

var cgroupFiles = []string{
	"cpu.max", "cpu.stat", "cpu.pressure", "memory.pressure", "io.pressure",
	"memory.current", "memory.max", "memory.peak", "memory.events", "memory.numa_stat",
	"cpuset.cpus", "cpuset.cpus.effective", "cpuset.cpus.exclusive.effective",
	"cpuset.cpus.partition", "cpuset.mems", "cpuset.mems.effective",
}

func read(path string, limit int) *string {
	raw, err := os.ReadFile(path)
	if err != nil {
		return nil
	}
	if len(raw) > limit {
		raw = raw[:limit]
	}
	s := string(raw)
	return &s
}

func trim(path string) *string {
	raw := read(path, 1000)
	if raw == nil {
		return nil
	}
	s := strings.TrimSpace(*raw)
	return &s
}

func integer(path string) *int {
	raw := trim(path)
	if raw == nil {
		return nil
	}
	v, err := strconv.Atoi(*raw)
	if err != nil {
		return nil
	}
	return &v
}

func cpuList(s string) []int {
	seen := map[int]bool{}
	for _, term := range strings.Split(strings.TrimSpace(s), ",") {
		if term == "" {
			continue
		}
		lo, hi := term, term
		if a, b, ok := strings.Cut(term, "-"); ok {
			lo, hi = a, b
		}
		start, e1 := strconv.Atoi(strings.TrimSpace(lo))
		end, e2 := strconv.Atoi(strings.TrimSpace(hi))
		if e1 != nil || e2 != nil || start < 0 || end < start || end > 4096 {
			continue
		}
		for cpu := start; cpu <= end; cpu++ {
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

func ownAffinity(pid int) ([]int, error) {
	var mask [128]byte
	_, _, code := syscall.RawSyscall(syscall.SYS_SCHED_GETAFFINITY,
		uintptr(pid), uintptr(len(mask)), uintptr(unsafe.Pointer(&mask[0])))
	if code != 0 {
		return nil, code
	}
	cpus := make([]int, 0)
	for i, bits := range mask {
		for bit := range 8 {
			if bits&(1<<bit) != 0 {
				cpus = append(cpus, i*8+bit)
			}
		}
	}
	return cpus, nil
}

func setAffinity(pid int, cpus []int) error {
	var mask [128]byte
	if len(cpus) == 0 {
		return fmt.Errorf("empty CPU affinity")
	}
	for _, cpu := range cpus {
		if cpu < 0 || cpu >= len(mask)*8 {
			return fmt.Errorf("CPU %d outside affinity bitset", cpu)
		}
		mask[cpu/8] |= 1 << (cpu % 8)
	}
	_, _, code := syscall.RawSyscall(syscall.SYS_SCHED_SETAFFINITY,
		uintptr(pid), uintptr(len(mask)), uintptr(unsafe.Pointer(&mask[0])))
	if code != 0 {
		return code
	}
	return nil
}

func readTopology(allowed []int) (map[int]topology.Place, map[int][]int) {
	nodes := map[int][]int{}
	nodeOf := map[int]int{}
	paths, _ := filepath.Glob("/sys/devices/system/node/node[0-9]*")
	for _, path := range paths {
		n, err := strconv.Atoi(strings.TrimPrefix(filepath.Base(path), "node"))
		if err != nil {
			continue
		}
		list := trim(filepath.Join(path, "cpulist"))
		if list != nil {
			nodes[n] = cpuList(*list)
			for _, cpu := range nodes[n] {
				nodeOf[cpu] = n
			}
		}
	}
	places := make(map[int]topology.Place, len(allowed))
	for _, cpu := range allowed {
		prefix := fmt.Sprintf("/sys/devices/system/cpu/cpu%d/topology", cpu)
		p := topology.Place{CPU: cpu,
			Package: integer(filepath.Join(prefix, "physical_package_id")),
			Core:    integer(filepath.Join(prefix, "core_id"))}
		if node, ok := nodeOf[cpu]; ok {
			p.Node = &node
		}
		if raw := trim(filepath.Join(prefix, "thread_siblings_list")); raw != nil {
			p.Siblings = cpuList(*raw)
		}
		places[cpu] = p
	}
	return places, nodes
}

type cpuTicks struct {
	User   uint64 `json:"user"`
	System uint64 `json:"system"`
	Idle   uint64 `json:"idle"`
	Steal  uint64 `json:"steal"`
}

type snapshot struct {
	UnixNS       int64              `json:"unix_ns"`
	CPUTicks     map[int]cpuTicks   `json:"cpu_ticks"`
	Interrupts   map[int]uint64     `json:"interrupts"`
	Pressure     map[string]*string `json:"pressure"`
	Cgroup       map[string]*string `json:"cgroup"`
	FrequencyKHz map[int]*string    `json:"frequency_khz"`
}

func makeSnapshot(cpus []int) snapshot {
	selected := map[int]bool{}
	for _, cpu := range cpus {
		selected[cpu] = true
	}
	ticks := map[int]cpuTicks{}
	if raw := read("/proc/stat", 1_000_000); raw != nil {
		for _, line := range strings.Split(*raw, "\n") {
			fields := strings.Fields(line)
			if len(fields) < 9 || !strings.HasPrefix(fields[0], "cpu") {
				continue
			}
			cpu, err := strconv.Atoi(strings.TrimPrefix(fields[0], "cpu"))
			if err != nil || !selected[cpu] {
				continue
			}
			v := make([]uint64, 8)
			valid := true
			for i := range v {
				v[i], err = strconv.ParseUint(fields[i+1], 10, 64)
				if err != nil {
					valid = false
					break
				}
			}
			if valid {
				ticks[cpu] = cpuTicks{User: v[0], System: v[2], Idle: v[3], Steal: v[7]}
			}
		}
	}
	interrupts := map[int]uint64{}
	if raw := read("/proc/interrupts", 1_000_000); raw != nil {
		lines := strings.Split(*raw, "\n")
		if len(lines) > 0 {
			columns := make([]int, 0)
			for _, field := range strings.Fields(lines[0]) {
				if cpu, err := strconv.Atoi(strings.TrimPrefix(field, "CPU")); err == nil {
					columns = append(columns, cpu)
				}
			}
			for _, line := range lines[1:] {
				fields := strings.Fields(line)
				if len(fields) < len(columns)+1 {
					continue
				}
				for i, cpu := range columns {
					if !selected[cpu] {
						continue
					}
					if value, err := strconv.ParseUint(fields[i+1], 10, 64); err == nil {
						interrupts[cpu] += value
					}
				}
			}
		}
	}
	cgroup := map[string]*string{}
	for _, name := range cgroupFiles {
		cgroup[name] = read(filepath.Join("/sys/fs/cgroup", name), 4096)
	}
	freq := map[int]*string{}
	for _, cpu := range cpus {
		freq[cpu] = trim(fmt.Sprintf(
			"/sys/devices/system/cpu/cpu%d/cpufreq/scaling_cur_freq", cpu))
	}
	return snapshot{UnixNS: time.Now().UnixNano(), CPUTicks: ticks,
		Interrupts: interrupts, Cgroup: cgroup,
		Pressure: map[string]*string{
			"cpu":    read("/proc/pressure/cpu", 4096),
			"memory": read("/proc/pressure/memory", 4096),
			"io":     read("/proc/pressure/io", 4096)}, FrequencyKHz: freq}
}

type irqRoute struct {
	SelectedCPUInterrupts uint64  `json:"selected_cpu_interrupts"`
	SMPAffinityList       *string `json:"smp_affinity_list"`
}

func irqRouting(cpus []int) map[string]irqRoute {
	raw := read("/proc/interrupts", 1_000_000)
	if raw == nil {
		return nil
	}
	lines := strings.Split(*raw, "\n")
	if len(lines) == 0 {
		return nil
	}
	columns := make([]int, 0)
	for _, field := range strings.Fields(lines[0]) {
		if cpu, err := strconv.Atoi(strings.TrimPrefix(field, "CPU")); err == nil {
			columns = append(columns, cpu)
		}
	}
	selected := map[int]bool{}
	for _, cpu := range cpus {
		selected[cpu] = true
	}
	type ranked struct {
		name  string
		count uint64
	}
	items := make([]ranked, 0)
	for _, line := range lines[1:] {
		fields := strings.Fields(line)
		if len(fields) < len(columns)+1 {
			continue
		}
		var count uint64
		for i, cpu := range columns {
			value, err := strconv.ParseUint(fields[i+1], 10, 64)
			if err != nil {
				count = 0
				break
			}
			if selected[cpu] {
				count += value
			}
		}
		if count > 0 {
			items = append(items, ranked{name: strings.TrimSuffix(fields[0], ":"), count: count})
		}
	}
	sort.Slice(items, func(i, j int) bool { return items[i].count > items[j].count })
	routes := map[string]irqRoute{}
	for _, item := range items[:min(len(items), 64)] {
		routes[item.name] = irqRoute{
			SelectedCPUInterrupts: item.count,
			SMPAffinityList:       trim(filepath.Join("/proc/irq", item.name, "smp_affinity_list")),
		}
	}
	return routes
}

type perfEvent struct {
	RawCount      uint64  `json:"raw_count"`
	TimeEnabledNS uint64  `json:"time_enabled_ns"`
	TimeRunningNS uint64  `json:"time_running_ns"`
	RunningFrac   float64 `json:"running_fraction"`
	ScaledCount   float64 `json:"scaled_count"`
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

func openPerf(pid int, config uint64) (int, *string) {
	var attr [64]byte
	binary.LittleEndian.PutUint32(attr[0:4], 0) // PERF_TYPE_HARDWARE
	binary.LittleEndian.PutUint32(attr[4:8], 64)
	binary.LittleEndian.PutUint64(attr[8:16], config)
	binary.LittleEndian.PutUint64(attr[32:40], 3) // enabled and running times
	binary.LittleEndian.PutUint64(attr[40:48], 2) // inherit child threads
	fd, _, code := syscall.RawSyscall6(syscall.SYS_PERF_EVENT_OPEN,
		uintptr(unsafe.Pointer(&attr[0])), uintptr(pid), ^uintptr(0), ^uintptr(0), 0, 0)
	if code != 0 {
		s := code.Error()
		return -1, &s
	}
	return int(fd), nil
}

func closePerf(fd int) (*perfEvent, *string) {
	if fd < 0 {
		return nil, nil
	}
	defer syscall.Close(fd)
	var raw [24]byte
	n, err := syscall.Read(fd, raw[:])
	if err != nil || n != len(raw) {
		s := fmt.Sprintf("perf_short_read:%d:%v", n, err)
		return nil, &s
	}
	count := binary.LittleEndian.Uint64(raw[0:8])
	enabled := binary.LittleEndian.Uint64(raw[8:16])
	running := binary.LittleEndian.Uint64(raw[16:24])
	if enabled == 0 || running > enabled {
		s := "perf_invalid_coverage"
		return nil, &s
	}
	fraction := float64(running) / float64(enabled)
	scaled := 0.0
	if running > 0 {
		scaled = float64(count) * float64(enabled) / float64(running)
	}
	return &perfEvent{RawCount: count, TimeEnabledNS: enabled,
		TimeRunningNS: running, RunningFrac: fraction, ScaledCount: scaled}, nil
}

func finishPerf(cyclesFD, instructionsFD int, cyclesErr, instructionsErr *string) perfPair {
	cycles, cReadErr := closePerf(cyclesFD)
	instructions, iReadErr := closePerf(instructionsFD)
	if cyclesErr == nil {
		cyclesErr = cReadErr
	}
	if instructionsErr == nil {
		instructionsErr = iReadErr
	}
	pair := perfPair{Cycles: cycles, Instructions: instructions,
		CyclesError: cyclesErr, InstructionsError: instructionsErr}
	if cycles != nil && instructions != nil && cycles.ScaledCount > 0 {
		v := instructions.ScaledCount / cycles.ScaledCount
		pair.IPCScaledEstimate = &v
		if cycles.RunningFrac >= .99 && instructions.RunningFrac >= .99 {
			pair.BothAtLeast99Pct = true
			pair.IPC = &v
		}
	}
	return pair
}

func processRSS(pid int) uint64 {
	raw := read(fmt.Sprintf("/proc/%d/status", pid), 50_000)
	if raw == nil {
		return 0
	}
	for _, line := range strings.Split(*raw, "\n") {
		fields := strings.Fields(line)
		if len(fields) >= 2 && fields[0] == "VmRSS:" {
			value, _ := strconv.ParseUint(fields[1], 10, 64)
			return value * 1024
		}
	}
	return 0
}

func threadMasks(pid int) (int, []string) {
	paths, _ := filepath.Glob(fmt.Sprintf("/proc/%d/task/*/status", pid))
	masks := make([]string, 0)
	for _, path := range paths {
		if raw := read(path, 50_000); raw != nil {
			for _, line := range strings.Split(*raw, "\n") {
				if strings.HasPrefix(line, "Cpus_allowed_list:") {
					masks = append(masks, strings.TrimSpace(strings.TrimPrefix(line, "Cpus_allowed_list:")))
					break
				}
			}
		}
	}
	sort.Strings(masks)
	return len(paths), slices.Compact(masks)
}

func numaPages(pid int) map[string]uint64 {
	pages := map[string]uint64{}
	raw := read(fmt.Sprintf("/proc/%d/numa_maps", pid), 1_000_000)
	if raw == nil {
		return pages
	}
	for _, field := range strings.Fields(*raw) {
		if len(field) < 3 || field[0] != 'N' {
			continue
		}
		node, value, ok := strings.Cut(field, "=")
		if !ok || len(node) < 2 {
			continue
		}
		if _, err := strconv.Atoi(node[1:]); err != nil {
			continue
		}
		if n, err := strconv.ParseUint(value, 10, 64); err == nil {
			pages[node] += n
		}
	}
	return pages
}
