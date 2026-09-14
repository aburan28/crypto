"""Exact block relaxation with nested prefix-independent column spans."""
import block
import reuse


class Basis(block.Linear):
    def __init__(self, f, parent=None):
        self.f, self.lifts, self.kernel = f, False, []
        self.rows = parent.rows.copy() if parent is not None else {}
        if parent is not None:
            f.word('cached_basis_field_word_copy', len(self.rows))
        self.ordered = sorted(self.rows.items(), reverse=True)

    def insert(self, columns):
        for value in columns:
            if len(self.rows) == self.f.m:
                break
            value, _ = self.reduce(value)
            if value:
                self.rows[value.bit_length() - 1] = value, 0
                self.ordered = sorted(self.rows.items(), reverse=True)


class Search(reuse.Search):
    def __init__(self, f, curve, v, target, deadline, batch_inverse=True):
        super().__init__(f, curve, v, target, deadline, batch_inverse)
        self.stats.update(cached_full_rank_skips=0, cached_basis_field_slots=0)

    def walk(self, circuit, b, h2, forbidden, free, prefix, rhs, columns):
        self.check()
        f, v = self.f, self.v
        if free == v.d:
            with self.phase('fixed_span_setup'):
                self.fixed = [Basis(f)]
                for i in range(v.d):
                    basis = Basis(f, self.fixed[-1])
                    basis.insert([circuit.zLinear[i]] + circuit.mixed[i])
                    self.fixed.append(basis)
                self.stats['cached_basis_field_slots'] = max(self.stats['cached_basis_field_slots'],
                                                              sum(len(x.rows) for x in self.fixed))
        if not free:
            yield from super().walk(circuit, b, h2, forbidden, free, prefix, rhs, columns)
            return
        with self.phase('block_rejection'):
            self.stats['block_tests'] += 1
            cached = self.fixed[free]
            if len(cached.rows) == f.m:
                self.stats['cached_full_rank_skips'] += 1
                self.stats['full_rank_blocks'] += 1
                mask = None
            else:
                matrix = Basis(f, cached)
                matrix.insert(columns)
                if len(matrix.rows) == f.m:
                    self.stats['full_rank_blocks'] += 1
                    mask = None
                else:
                    mask = matrix.separator(rhs)
            if mask is not None:
                removed = (1 << free) - sum(x >> free == prefix >> free for x in forbidden)
                self.stats['pruned_admissible_branches'] += removed
                self.stats['rejected_blocks'] += 1
                bucket = self.stats['rejections_by_free_bits'].setdefault(str(free), {'blocks': 0, 'admissible_branches': 0})
                bucket['blocks'] += 1
                bucket['admissible_branches'] += removed
                self.certificates.append({'b': f.toCoords(b), 'h2': f.toCoords(h2),
                                          'prefix': prefix, 'free_bits': free, 'separator': mask})
                return
        bit = free - 1
        yield from self.walk(circuit, b, h2, forbidden, bit, prefix, rhs, columns)
        with self.phase('partition_updates'):
            next_columns = [f.add(a, delta) for a, delta in zip(columns, circuit.mixed[bit])]
            next_rhs = f.add(rhs, circuit.zLinear[bit])
        yield from self.walk(circuit, b, h2, forbidden, bit, prefix | (1 << bit), next_rhs, next_columns)
