package topology

import (
	"slices"
	"sort"
	"strconv"
)

// Place is the topology visible to the guest for one allowed logical CPU.
// Missing package/core IDs remain unknown and never count as proved cores.
type Place struct {
	CPU      int   `json:"cpu"`
	Package  *int  `json:"package"`
	Core     *int  `json:"core"`
	Node     *int  `json:"node"`
	Siblings []int `json:"siblings"`
}

type Selection struct {
	Allowed              []int  `json:"allowed"`
	CPUs                 []int  `json:"cpus"`
	KnownPhysicalCores   int    `json:"known_physical_cores"`
	CoreIdentityComplete bool   `json:"core_identity_complete"`
	WholeSMT             bool   `json:"whole_smt"`
	SingleNUMANode       bool   `json:"single_numa_node"`
	Mode                 string `json:"mode"`
}

type group struct {
	cpus []int
	node *int
}

// Select prefers complete SMT groups on one NUMA node, then uses at most the
// visible siblings of physicalLimit observed cores. It never assumes that a
// Modal physical-core request exposes both SMT siblings in sched_getaffinity.
func Select(allowed []int, places map[int]Place, physicalLimit, logicalLimit int) Selection {
	allowed = slices.Clone(allowed)
	sort.Ints(allowed)
	allowed = slices.Compact(allowed)
	result := Selection{Allowed: allowed, Mode: "no_eligible_cpus"}
	if len(allowed) == 0 || physicalLimit <= 0 || logicalLimit <= 0 {
		return result
	}
	allowedSet := make(map[int]bool, len(allowed))
	for _, cpu := range allowed {
		allowedSet[cpu] = true
	}
	groups := make([]group, 0)
	index := make(map[string]int)
	known := true
	for _, cpu := range allowed {
		p, ok := places[cpu]
		key := "unknown:" + strconv.Itoa(cpu)
		if ok && p.Package != nil && p.Core != nil {
			key = "core:" + strconv.Itoa(*p.Package) + ":" + strconv.Itoa(*p.Core)
		} else {
			known = false
		}
		if i, exists := index[key]; exists {
			groups[i].cpus = append(groups[i].cpus, cpu)
		} else {
			index[key] = len(groups)
			groups = append(groups, group{cpus: []int{cpu}, node: p.Node})
		}
	}

	complete := func(g group) bool {
		members := make(map[int]bool, len(g.cpus))
		for _, cpu := range g.cpus {
			members[cpu] = true
		}
		for _, cpu := range g.cpus {
			siblings := places[cpu].Siblings
			if len(siblings) == 0 || !slices.Contains(siblings, cpu) {
				return false
			}
			for _, sibling := range siblings {
				if !members[sibling] {
					return false
				}
			}
		}
		return true
	}
	filter := func(in []group, wholeOnly bool) []group {
		out := make([]group, 0, len(in))
		for _, g := range in {
			if !wholeOnly || complete(g) {
				out = append(out, g)
			}
		}
		return out
	}
	nodeGroups := func(node int) []group {
		out := make([]group, 0)
		for _, g := range groups {
			if g.node != nil && *g.node == node {
				out = append(out, g)
			}
		}
		return out
	}
	// Prefer a single node of complete core groups, then complete groups
	// across nodes. A partial affinity mask remains usable as a diagnostic.
	var chosen []group
	nodeIDs := make([]int, 0)
	for _, g := range groups {
		if g.node != nil && !slices.Contains(nodeIDs, *g.node) {
			nodeIDs = append(nodeIDs, *g.node)
		}
	}
	sort.Ints(nodeIDs)
	for _, node := range nodeIDs {
		if completeOnNode := filter(nodeGroups(node), true); len(completeOnNode) >= physicalLimit {
			chosen = completeOnNode
			break
		}
	}
	if chosen == nil {
		if completeAll := filter(groups, true); len(completeAll) >= physicalLimit {
			chosen = completeAll
		}
	}
	if chosen == nil {
		for _, node := range nodeIDs {
			if onNode := nodeGroups(node); len(onNode) >= physicalLimit {
				chosen = onNode
				break
			}
		}
	}
	if chosen == nil {
		chosen = groups
	}

	selectedGroups := chosen[:min(len(chosen), physicalLimit)]
	for _, g := range selectedGroups {
		for _, cpu := range g.cpus[:min(len(g.cpus), 2)] {
			if len(result.CPUs) < logicalLimit {
				result.CPUs = append(result.CPUs, cpu)
			}
		}
	}
	result.CoreIdentityComplete = known
	if known {
		result.KnownPhysicalCores = len(selectedGroups)
	}
	selectedSet := make(map[int]bool, len(result.CPUs))
	for _, cpu := range result.CPUs {
		selectedSet[cpu] = true
	}
	result.WholeSMT = len(result.CPUs) > 0
	for _, cpu := range result.CPUs {
		p := places[cpu]
		if len(p.Siblings) == 0 || !slices.Contains(p.Siblings, cpu) {
			result.WholeSMT = false
			break
		}
		for _, sibling := range p.Siblings {
			if !allowedSet[sibling] || !selectedSet[sibling] {
				result.WholeSMT = false
				break
			}
		}
	}
	result.SingleNUMANode = len(selectedGroups) > 0 && selectedGroups[0].node != nil
	if result.SingleNUMANode {
		for _, g := range selectedGroups[1:] {
			if g.node == nil || *g.node != *selectedGroups[0].node {
				result.SingleNUMANode = false
				break
			}
		}
	}
	result.Mode = "observed_affinity_cores"
	if result.WholeSMT && result.KnownPhysicalCores == physicalLimit &&
		len(result.CPUs) == logicalLimit {
		result.Mode = "complete_smt_groups"
	}
	return result
}
