import unittest
from oracle import InvalidEvidence
from test_tournament import rows
from tournament import comparison


class PolicyComparisonTests(unittest.TestCase):
    def data(self):
        data = rows(.6)
        for row in data:
            if row['arm'] == 'candidate':
                row['certificate']['factor_base_sha256'] = 'different-declared-support'
        return data

    def test_support_change_requires_explicit_policy_mode(self):
        with self.assertRaises(InvalidEvidence):
            comparison(self.data(), 'candidate', draws=100)
        self.assertTrue(comparison(self.data(), 'candidate', draws=100, match_support=False)['eligible'])

    def test_policy_still_rejects_changed_support_within_an_arm(self):
        data = self.data()
        next(row for row in data if row['arm'] == 'candidate')['certificate']['factor_base_sha256'] = 'unstable'
        with self.assertRaises(InvalidEvidence):
            comparison(data, 'candidate', draws=100, match_support=False)

    def test_policy_still_rejects_changed_public_targets(self):
        data = self.data()
        next(row for row in data if row['arm'] == 'candidate')['case_sha256'] = 'different-target'
        with self.assertRaises(InvalidEvidence):
            comparison(data, 'candidate', draws=100, match_support=False)


if __name__ == '__main__':
    unittest.main()
