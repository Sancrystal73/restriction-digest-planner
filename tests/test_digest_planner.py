"""Five test cases for the Restriction Digest Planner.

Written with unittest, the testing framework from class. Run from
the repository root with:

    python -m unittest discover -s tests
"""

import os
import unittest

from digest_planner.compare import distinguishable, separation_score
from digest_planner.digest import (
    cut_positions,
    find_sites,
    fragment_lengths,
    insertion_sort,
)
from digest_planner.enzymes import (
    Enzyme,
    is_palindromic,
    load_enzymes,
    reverse_complement,
)
from digest_planner.planner import DigestPlanner
from digest_planner.sequence_io import Plasmid

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "data")


class TestDigestPlanner(unittest.TestCase):

    def test_1_reverse_complement_and_palindrome(self):
        """Reverse complements and the palindrome check."""
        self.assertEqual(reverse_complement("GAATTC"), "GAATTC")
        self.assertEqual(reverse_complement("GGTCTC"), "GAGACC")
        self.assertTrue(is_palindromic("GAATTC"))
        self.assertFalse(is_palindromic("GGTCTC"))
        # Caret parsing: the cut offset is the position of '^'.
        pst1 = Enzyme.from_caret("PstI", "CTGCA^G")
        self.assertEqual(pst1.site, "CTGCAG")
        self.assertEqual(pst1.cut_offset, 5)

    def test_2_find_sites_and_sort(self):
        """Site search (linear and origin-spanning) and insertion sort."""
        self.assertEqual(find_sites("AAGAATTCAAGAATTC", "GAATTC", False),
                         [2, 10])
        # The EcoRI site starts at position 10 and wraps the origin.
        self.assertEqual(find_sites("ATTCAAAAAAGA", "GAATTC", True), [10])
        # EcoRI cuts one base in: (10 + 1) mod 12 = 11.
        plasmid = Plasmid("wrap", "ATTCAAAAAAGA", True)
        ecori = Enzyme("EcoRI", "GAATTC", 1)
        self.assertEqual(cut_positions(plasmid, [ecori]), [11])
        self.assertEqual(insertion_sort([2400, 15, 900, 15]),
                         [15, 15, 900, 2400])
        self.assertEqual(insertion_sort([]), [])
        self.assertEqual(insertion_sort([7]), [7])

    def test_3_fragment_lengths(self):
        """Fragment lengths from cut positions, including edge cases."""
        # Two cuts: the gaps, plus the fragment wrapping the origin.
        self.assertEqual(fragment_lengths([100, 900], 3000), [800, 2200])
        # Zero cuts: uncut. One cut: one fragment of full length.
        self.assertEqual(fragment_lengths([], 3000), [])
        self.assertEqual(fragment_lengths([1234], 3000), [3000])
        # Fragments always add up to the plasmid length.
        self.assertEqual(sum(fragment_lengths([100, 900, 2500], 3000)), 3000)

    def test_4_version_1_vs_version_2(self):
        """The threshold improvement: near-identical bands are not a hit."""
        a = [1000, 2000]
        b = [1020, 1980]
        # Version 1 counts any difference at all.
        self.assertTrue(distinguishable(a, b, None))
        # Version 2 at 0.05 does not: the score is only 0.0196.
        self.assertAlmostEqual(separation_score(a, b), 0.0196, places=4)
        self.assertFalse(distinguishable(a, b, 0.05))
        # A clearly different pair still counts at 0.10.
        self.assertAlmostEqual(separation_score([800, 2200], [1400, 1600]),
                               0.43, places=2)
        self.assertTrue(distinguishable([800, 2200], [1400, 1600], 0.10))

    def test_5_synthetic_pair_end_to_end(self):
        """The full program on the synthetic pair (worked example)."""
        enzymes = load_enzymes(os.path.join(DATA, "enzymes.csv"))
        plasmid_a = Plasmid.from_fasta(os.path.join(DATA, "plasmid_A.fasta"))
        plasmid_b = Plasmid.from_fasta(os.path.join(DATA, "plasmid_B.fasta"))
        planner = DigestPlanner(plasmid_a, plasmid_b, enzymes, 0.10)
        planner.run()

        # Single digests only open the circle: not distinguishing.
        single = planner.evaluate(("EcoRI",))
        self.assertEqual(single.frags_a, [3000])
        self.assertEqual(single.frags_b, [3000])
        self.assertFalse(single.distinguishable)

        # The double digest distinguishes, with the expected fragments.
        double = planner.evaluate(("BamHI", "EcoRI"))
        self.assertEqual(sorted(double.frags_a), [800, 2200])
        self.assertEqual(sorted(double.frags_b), [1400, 1600])
        self.assertTrue(double.distinguishable)
        self.assertAlmostEqual(double.score, 0.43, places=2)

        # It is the only distinguishing digest out of all 78.
        best = planner.top(5)
        self.assertEqual(len(best), 1)
        self.assertEqual(best[0].scheme, ("BamHI", "EcoRI"))


if __name__ == "__main__":
    unittest.main()
