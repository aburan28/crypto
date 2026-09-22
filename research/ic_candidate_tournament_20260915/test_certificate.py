"""The index-calculus admission requirement: every target's logarithm must
carry the descent relation it was derived from, and the checker must reject
anything else. The report below is a real output of the round-0006 baseline
worker (degree 13, two public-hash targets, no planted scalar)."""
import copy
import unittest

from oracle import Curve, InvalidEvidence, verify

FIXTURE = {"cofactor": "4", "curve_a": 0, "degree": 13, "generator": ["4793", "7108"], "group_order": "8012", "irreducible": {"degree": 13, "low_terms": [0, 1, 3, 4]}, "lambda": "89", "subgroup_order": "2003", "target_scalar_constructed": False, "target_seeds": [1, 2], "targets": [["6260", "4377"], ["5985", "7591"]]}

REPORT = {"accepted_relations": 8, "column_logs": [{"log": "1841", "point": ["22", "7783"]}, {"log": "1573", "point": ["48", "3480"]}, {"log": "763", "point": ["118", "3210"]}, {"log": "104", "point": ["521", "1050"]}, {"log": "1724", "point": ["755", "7558"]}, {"log": "301", "point": ["931", "2283"]}, {"log": "1897", "point": ["949", "3443"]}], "columns": 7, "duplicate_relations": 0, "factor_base": [["44", "564"], ["44", "536"], ["50", "6388"], ["50", "6342"], ["74", "5509"], ["74", "5583"], ["130", "4488"], ["130", "4362"], ["154", "4606"], ["154", "4452"], ["156", "7223"], ["156", "7339"], ["266", "3878"], ["266", "3628"], ["358", "1748"], ["358", "1970"], ["370", "3552"], ["370", "3218"], ["655", "516"], ["655", "139"], ["771", "880"], ["771", "115"], ["913", "6124"], ["913", "5245"], ["957", "5816"], ["957", "5381"], ["1011", "7532"], ["1011", "7839"], ["1013", "4077"], ["1013", "3096"], ["1015", "6607"], ["1015", "6712"], ["1104", "1648"], ["1104", "544"], ["1124", "3378"], ["1124", "2390"], ["1156", "7059"], ["1156", "7959"], ["1284", "7015"], ["1284", "7779"], ["1374", "7561"], ["1374", "6359"], ["1412", "6415"], ["1412", "7307"], ["1468", "4187"], ["1468", "5607"], ["1557", "1073"], ["1557", "548"], ["1605", "4659"], ["1605", "5238"], ["1631", "6535"], ["1631", "8152"], ["1741", "2789"], ["1741", "3112"], ["1759", "5068"], ["1759", "5395"], ["1813", "5584"], ["1813", "4805"], ["2001", "2147"], ["2001", "4018"], ["2110", "2335"], ["2110", "289"], ["2122", "4371"], ["2122", "6489"], ["2544", "3304"], ["2544", "1304"], ["2819", "6911"], ["2819", "4604"], ["2999", "5952"], ["2999", "7415"], ["3013", "4127"], ["3013", "7130"], ["3029", "1787"], ["3029", "3374"], ["3098", "5429"], ["3098", "6447"], ["3104", "3724"], ["3104", "684"], ["3348", "5676"], ["3348", "6968"], ["3400", "6636"], ["3400", "5284"], ["3454", "3788"], ["3454", "946"], ["3494", "3274"], ["3494", "364"], ["3607", "2205"], ["3607", "1674"], ["3881", "1076"], ["3881", "2845"], ["4015", "4276"], ["4015", "7963"], ["4081", "2177"], ["4081", "1904"], ["4164", "5413"], ["4164", "1377"], ["4200", "5505"], ["4200", "1513"], ["4224", "1634"], ["4224", "5858"], ["4384", "7749"], ["4384", "3941"], ["4570", "6933"], ["4570", "2767"], ["4943", "6038"], ["4943", "1241"], ["5057", "5705"], ["5057", "1416"], ["5061", "8140"], ["5061", "3081"], ["5152", "6270"], ["5152", "3166"], ["5240", "468"], ["5240", "5548"], ["5324", "8134"], ["5324", "2826"], ["5526", "1524"], ["5526", "4194"], ["5542", "6210"], ["5542", "3556"], ["5596", "3957"], ["5596", "6825"], ["5771", "5961"], ["5771", "450"], ["5787", "7930"], ["5787", "2145"], ["5791", "3140"], ["5791", "6875"], ["5833", "1571"], ["5833", "4330"], ["6085", "5933"], ["6085", "232"], ["6089", "5587"], ["6089", "538"], ["6105", "2955"], ["6105", "7250"], ["6226", "1724"], ["6226", "7918"], ["6252", "6884"], ["6252", "648"], ["6544", "7751"], ["6544", "2007"], ["6638", "6012"], ["6638", "3730"], ["6673", "6407"], ["6673", "790"], ["6739", "4842"], ["6739", "2233"], ["6879", "4115"], ["6879", "2764"], ["6907", "818"], ["6907", "6601"], ["7180", "878"], ["7180", "8034"], ["7204", "3912"], ["7204", "4972"], ["7232", "7548"], ["7232", "316"], ["7234", "6495"], ["7234", "1309"], ["7266", "2192"], ["7266", "5362"], ["7296", "7136"], ["7296", "1888"], ["7298", "4272"], ["7298", "3122"], ["7815", "348"], ["7815", "8155"], ["7921", "7743"], ["7921", "206"], ["7951", "542"], ["7951", "7441"], ["8005", "2950"], ["8005", "5315"], ["8033", "3049"], ["8033", "5256"], ["8067", "2948"], ["8067", "5127"], ["8091", "6301"], ["8091", "1798"], ["8117", "3770"], ["8117", "4367"]], "fixture": {"cofactor": "4", "curve_a": 0, "degree": 13, "generator": ["4793", "7108"], "group_order": "8012", "irreducible": {"degree": 13, "low_terms": [0, 1, 3, 4]}, "lambda": "89", "subgroup_order": "2003", "target_scalar_constructed": False, "target_seeds": [1, 2], "targets": [["6260", "4377"], ["5985", "7591"]]}, "mode": "ic", "rejected_relations": 0, "relations": [{"a": 1998, "points": [0, 2, 38], "trial": 0}, {"a": 641, "points": [0, 14, 48], "trial": 1}, {"a": 856, "points": [0, 10, 135], "trial": 2}, {"a": 1026, "points": [0, 9, 166], "trial": 3}, {"a": 335, "points": [0, 25, 174], "trial": 4}, {"a": 12, "points": [0, 7, 36], "trial": 5}, {"a": 152, "points": [0, 2, 109], "trial": 6}, {"a": 1117, "points": [0, 14, 175], "trial": 7}], "schema_version": 1, "solutions": [{"index": 0, "recovered": "1621", "relation": {"a": 12, "b": 128, "points": [0, 9, 25]}, "trials": 1}, {"index": 1, "recovered": "1301", "relation": {"a": 12, "b": 128, "points": [0, 2, 6]}, "trials": 1}], "solve_attempts": 2, "sparse_report": {"attempts": 0, "core_dimension": 0, "core_nonzeros": 0, "filter": {"columns_in": 7, "columns_out": 0, "dependent_rows_dropped": 1, "duplicates_removed": 0, "excess_rows_removed": 0, "merged_columns": 5, "nonzeros_in": 23, "nonzeros_out": 0, "rows_in": 8, "rows_out": 0, "singletons_removed": 2, "uncovered_columns": 0}, "reconstructed_columns": 7, "wiedemann": None}, "status": "complete", "summands": 3, "trials": 8}


class DescentCertificateTests(unittest.TestCase):
    def test_certified_report_verifies_every_target(self):
        proof = verify(copy.deepcopy(REPORT), FIXTURE, expected_mode='ic', summands=3)
        self.assertEqual(proof['certified_descents'], 2)
        self.assertEqual(proof['degenerate_descents'], 0)
        self.assertEqual(proof['solutions'], ['1621', '1301'])

    def rejects(self, mutate, message):
        report = copy.deepcopy(REPORT)
        mutate(report)
        with self.assertRaises(InvalidEvidence, msg=message) as caught:
            verify(report, FIXTURE, expected_mode='ic', summands=3)
        self.assertIn(message, str(caught.exception))

    def test_logarithm_without_relation_is_not_index_calculus(self):
        self.rejects(lambda r: r['solutions'][0].pop('relation'), 'not certified as index calculus')

    def test_relation_must_hold_in_the_group(self):
        self.rejects(lambda r: r['solutions'][0]['relation'].update(a=13), 'does not hold in the group')
        self.rejects(lambda r: r['solutions'][0]['relation'].update(points=[0, 9, 26]), 'does not hold in the group')
        self.rejects(lambda r: r['solutions'][0]['relation'].update(points=[]), 'does not hold in the group')

    def test_relation_belongs_to_its_own_target(self):
        def swap(r):
            first, second = r['solutions']
            first['relation'], second['relation'] = second['relation'], first['relation']
        self.rejects(swap, 'does not hold in the group')

    def test_scalar_and_probe_coefficients_are_checked(self):
        self.rejects(lambda r: r['solutions'][0]['relation'].update(b=0), 'invalid descent scalars')
        self.rejects(lambda r: r['solutions'][0].update(recovered='1622'), 'incorrect scalar')

    def test_equivalent_coefficients_for_the_same_probe_remain_certified(self):
        # [a]G + [b]Q = [a']G + [b']Q whenever a + b·d = a' + b'·d: the same probe
        # point decomposed over the same base points is the same index-calculus
        # derivation, so the certificate still verifies.
        report = copy.deepcopy(REPORT)
        rel = report['solutions'][0]['relation']
        rel['a'], rel['b'] = (rel['a'] + 1621) % 2003, (rel['b'] + 2002) % 2003
        proof = verify(report, FIXTURE, expected_mode='ic', summands=3)
        self.assertEqual(proof['certified_descents'], 2)

    def test_degenerate_relation_is_only_accepted_when_the_probe_is_infinity(self):
        self.rejects(lambda r: r['solutions'][1]['relation'].update(points=[]), 'does not hold in the group')


class DegreeBoundTests(unittest.TestCase):
    """The checker's degree bound mirrors the collector's MAX_DEGREE.

    It is not a width limit -- the arithmetic here is Python integers over the
    fixture's own irreducible polynomial. It exists so the checker refuses
    exactly what `koblitz_tiny_ic` refuses. Round 0021 widened the pair table's
    packed coordinates from `u32` to `u64`, which moved that ceiling from 31 to
    61, and this bound followed.

    The bound is asserted by WHICH failure a degree produces, not by whether
    one does. `Curve.__init__` goes on to check the group order against the
    Frobenius trace, the cofactor, the generator and the eigenvalue, so a
    fixture with only its degree swapped fails on the group order no matter
    what the bound says. A degree the bound rejects therefore has to fail with
    'unsupported degree'; a degree it accepts has to fail with something else,
    which is what proves it got past the bound.
    """

    def failure_for(self, degree):
        fixture = dict(FIXTURE, degree=degree,
                       irreducible={'degree': degree, 'low_terms': [0, 1, 3, 4]})
        try:
            Curve(fixture)
        except InvalidEvidence as exc:
            return str(exc)
        return ''

    def test_accepts_the_degrees_the_collector_now_runs(self):
        for degree in (5, 13, 31, 37, 41, 61):
            self.assertNotIn('unsupported degree', self.failure_for(degree), f'degree {degree}')

    def test_refuses_above_the_collector_ceiling(self):
        for degree in (63, 65, 127):
            self.assertIn('unsupported degree', self.failure_for(degree), f'degree {degree}')

    def test_still_refuses_even_degrees_and_tiny_fields(self):
        for degree in (3, 4, 12, 36):
            self.assertIn('unsupported degree', self.failure_for(degree), f'degree {degree}')

    def test_the_round_0020_fixture_still_verifies_unchanged(self):
        """The amendment is additive: a degree-13 report reads exactly as before."""
        self.assertEqual(verify(copy.deepcopy(REPORT), FIXTURE, expected_mode='ic',
                                summands=3)['solutions'], ['1621', '1301'])

if __name__ == '__main__':
    unittest.main()
