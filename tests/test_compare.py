"""Tests for compare.py: T11, T12."""

import unittest

from digest_planner.compare import (
    differs_exact,
    distinguishable,
    separation_score,
)


class TestVersion1(unittest.TestCase):
    def test_t11_distinguishing(self):
        self.assertTrue(differs_exact([800, 2200], [1400, 1600]))

    def test_t11_not_distinguishing(self):
        self.assertFalse(differs_exact([1500, 1500], [1500, 1500]))

    def test_order_does_not_matter(self):
        self.assertFalse(differs_exact([800, 2200], [2200, 800]))


class TestVersion2(unittest.TestCase):
    def test_t12_threshold(self):
        a = [1000, 2000]
        b = [1020, 1980]
        # Version 1: any difference counts.
        self.assertTrue(distinguishable(a, b, threshold=None))
        # Version 2 at 0.05: score is 0.0196, below threshold.
        self.assertAlmostEqual(separation_score(a, b), 0.0196, places=4)
        self.assertFalse(distinguishable(a, b, threshold=0.05))

    def test_t13_score(self):
        score = separation_score([800, 2200], [1400, 1600])
        self.assertAlmostEqual(score, 0.4286, places=4)
        self.assertTrue(distinguishable([800, 2200], [1400, 1600], threshold=0.10))

    def test_identical_scores_zero(self):
        self.assertEqual(separation_score([1500, 1500], [1500, 1500]), 0.0)

    def test_one_uncut(self):
        self.assertTrue(distinguishable([], [3000], threshold=0.10))
        self.assertTrue(distinguishable([], [3000], threshold=None))
        self.assertEqual(separation_score([], [3000]), 0.0)
        self.assertFalse(distinguishable([], [], threshold=0.10))


if __name__ == "__main__":
    unittest.main()
