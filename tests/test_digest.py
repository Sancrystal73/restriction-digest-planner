"""Tests for digest.py: T4-T10, plus input validation T15."""

import random
import tempfile
import unittest
from pathlib import Path

from digest_planner.digest import (
    cut_positions,
    find_sites,
    fragment_lengths,
    selection_sort,
)
from digest_planner.enzymes import Enzyme
from digest_planner.sequence_io import Plasmid


class TestFindSites(unittest.TestCase):
    def test_t4_linear(self):
        self.assertEqual(
            find_sites("AAGAATTCAAGAATTC", "GAATTC", circular=False), [2, 10]
        )

    def test_t5_spanning_origin(self):
        self.assertEqual(
            find_sites("ATTCAAAAAAGA", "GAATTC", circular=True), [10]
        )
        # EcoRI cut at (10 + 1) mod 12 = 11.
        plasmid = Plasmid("wrap", "ATTCAAAAAAGA", circular=True)
        ecori = Enzyme("EcoRI", "GAATTC", 1)
        self.assertEqual(cut_positions(plasmid, [ecori]), [11])

    def test_no_site(self):
        self.assertEqual(find_sites("AAAAAA", "GAATTC", circular=False), [])
        self.assertEqual(find_sites("AAAAAA", "GAATTC", circular=True), [])


class TestSelectionSort(unittest.TestCase):
    def test_t6_typical(self):
        self.assertEqual(selection_sort([2400, 15, 900, 15]), [15, 15, 900, 2400])

    def test_t6_empty_and_single(self):
        self.assertEqual(selection_sort([]), [])
        self.assertEqual(selection_sort([7]), [7])

    def test_does_not_mutate_input(self):
        values = [3, 1, 2]
        self.assertEqual(selection_sort(values), [1, 2, 3])
        self.assertEqual(values, [3, 1, 2])


class TestFragmentLengths(unittest.TestCase):
    def test_t7_circular_two_cuts(self):
        self.assertEqual(fragment_lengths([100, 900], 3000), [800, 2200])

    def test_t8_zero_and_one_cut(self):
        self.assertEqual(fragment_lengths([], 3000), [])
        self.assertEqual(fragment_lengths([1234], 3000), [3000])

    def test_fragments_sum_to_length(self):
        self.assertEqual(sum(fragment_lengths([100, 900, 2500], 3000)), 3000)

    def test_unsorted_input(self):
        self.assertEqual(fragment_lengths([900, 100], 3000), [800, 2200])


class TestCutPositions(unittest.TestCase):
    def test_t9_duplicate_cut_removed(self):
        # Two enzymes with the same site and cut position.
        plasmid = Plasmid("dup", "AAGAATTC" + "A" * 92, circular=True)
        e1 = Enzyme("EcoRI", "GAATTC", 1)
        e2 = Enzyme("EcoRI-copy", "GAATTC", 1)
        cuts = cut_positions(plasmid, [e1, e2])
        self.assertEqual(cuts, [3])
        self.assertEqual(fragment_lengths(cuts, plasmid.length), [100])

    def test_t10_fragment_sum_invariant(self):
        rng = random.Random(123)
        enzymes = [
            Enzyme("EcoRI", "GAATTC", 1),
            Enzyme("BamHI", "GGATCC", 1),
            Enzyme("HindIII", "AAGCTT", 1),
        ]
        for _ in range(100):
            seq = "".join(rng.choice("ACGT") for _ in range(3000))
            plasmid = Plasmid("rand", seq, circular=True)
            for enzyme in enzymes:
                cuts = cut_positions(plasmid, [enzyme])
                frags = fragment_lengths(cuts, plasmid.length)
                if cuts:
                    self.assertEqual(sum(frags), 3000)
            cuts = cut_positions(plasmid, enzymes)
            frags = fragment_lengths(cuts, plasmid.length)
            if cuts:
                self.assertEqual(sum(frags), 3000)


class TestInputValidation(unittest.TestCase):
    def test_t15_invalid_base_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.fasta"
            path.write_text(">bad\nACGTNACGT\n")
            with self.assertRaises(ValueError) as ctx:
                Plasmid.from_fasta(path)
            self.assertIn("N", str(ctx.exception))
            self.assertIn("4", str(ctx.exception))  # position of N

    def test_lowercase_accepted(self):
        plasmid = Plasmid("low", "gaattc", circular=False)
        self.assertEqual(plasmid.sequence, "GAATTC")


if __name__ == "__main__":
    unittest.main()
