package topology

import (
	"reflect"
	"testing"
)

func pointer(n int) *int { return &n }

func interleaved(cores int, visibleSiblings bool) ([]int, map[int]Place) {
	allowed := make([]int, 0, cores*2)
	places := make(map[int]Place)
	for core := range cores {
		for thread := range 2 {
			if thread == 1 && !visibleSiblings {
				continue
			}
			cpu := core + thread*cores
			allowed = append(allowed, cpu)
			places[cpu] = Place{CPU: cpu, Package: pointer(0), Core: pointer(core),
				Node: pointer(0), Siblings: []int{core, core + cores}}
		}
	}
	return allowed, places
}

func TestCompleteInterleavedSMT(t *testing.T) {
	allowed, places := interleaved(8, true)
	got := Select(allowed, places, 8, 16)
	if len(got.CPUs) != 16 || got.KnownPhysicalCores != 8 || !got.WholeSMT ||
		!got.SingleNUMANode || got.Mode != "complete_smt_groups" {
		t.Fatalf("complete SMT selection: %+v", got)
	}
	if !reflect.DeepEqual(got.CPUs[:4], []int{0, 8, 1, 9}) {
		t.Fatalf("unexpected interleaved selection: %v", got.CPUs)
	}
}

func TestPartialSMTAffinityKeepsObservedCores(t *testing.T) {
	allowed, places := interleaved(8, false)
	got := Select(allowed, places, 8, 16)
	if !reflect.DeepEqual(got.CPUs, []int{0, 1, 2, 3, 4, 5, 6, 7}) ||
		got.KnownPhysicalCores != 8 || got.WholeSMT ||
		got.Mode != "observed_affinity_cores" {
		t.Fatalf("partial SMT selection: %+v", got)
	}
}

func TestMoreGuestCoresDoNotIncreaseEnvelope(t *testing.T) {
	allowed, places := interleaved(16, true)
	got := Select(allowed, places, 8, 16)
	if len(got.CPUs) != 16 || got.KnownPhysicalCores != 8 || !got.WholeSMT {
		t.Fatalf("oversized topology selection: %+v", got)
	}
}

func TestMissingCoreIDsRemainUnverified(t *testing.T) {
	allowed := make([]int, 12)
	places := make(map[int]Place)
	for cpu := range allowed {
		allowed[cpu] = cpu
		places[cpu] = Place{CPU: cpu, Siblings: []int{cpu}}
	}
	got := Select(allowed, places, 8, 16)
	if len(got.CPUs) != 8 || got.KnownPhysicalCores != 0 || got.CoreIdentityComplete {
		t.Fatalf("missing core IDs incorrectly verified: %+v", got)
	}
}

func TestNoAffinityEmitsNoEligibleSelection(t *testing.T) {
	got := Select(nil, nil, 8, 16)
	if len(got.CPUs) != 0 || got.Mode != "no_eligible_cpus" {
		t.Fatalf("empty selection: %+v", got)
	}
}
