"""Tests for enzymes.py: T1, T2, T3."""

import csv
import tempfile
import unittest
import warnings
from pathlib import Path

from digest_planner.enzymes import (
    Enzyme,
    is_palindromic,
    load_enzymes,
    reverse_complement,
)


class TestReverseComplement(unittest.TestCase):
    def test_t1_palindromic_site(self):
        self.assertEqual(reverse_complement("GAATTC"), "GAATTC")

    def test_t1_non_palindromic_site(self):
        self.assertEqual(reverse_complement("GGTCTC"), "GAGACC")

    def test_is_palindromic(self):
        self.assertTrue(is_palindromic("GAATTC"))
        self.assertFalse(is_palindromic("GGTCTC"))


class TestCaretParsing(unittest.TestCase):
    def test_t3_psti(self):
        enzyme = Enzyme.from_caret("PstI", "CTGCA^G")
        self.assertEqual(enzyme.site, "CTGCAG")
        self.assertEqual(enzyme.cut_offset, 5)

    def test_t3_ecori(self):
        enzyme = Enzyme.from_caret("EcoRI", "G^AATTC")
        self.assertEqual(enzyme.site, "GAATTC")
        self.assertEqual(enzyme.cut_offset, 1)

    def test_missing_caret_raises(self):
        with self.assertRaises(ValueError):
            Enzyme.from_caret("Bad", "GAATTC")


class TestLoader(unittest.TestCase):
    def test_t2_non_palindromic_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "enzymes.csv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["name", "site"])
                writer.writerow(["EcoRI", "G^AATTC"])
                writer.writerow(["BsaI", "GGTCT^C"])  # non-palindromic GGTCTC
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                enzymes = load_enzymes(path)
            self.assertIn("EcoRI", enzymes)
            self.assertNotIn("BsaI", enzymes)
            self.assertTrue(any("BsaI" in str(w.message) for w in caught))

    def test_load_real_table(self):
        root = Path(__file__).resolve().parent.parent
        enzymes = load_enzymes(root / "data" / "enzymes.csv")
        self.assertEqual(len(enzymes), 12)
        self.assertEqual(enzymes["EcoRI"].site, "GAATTC")
        self.assertEqual(enzymes["EcoRI"].cut_offset, 1)
        self.assertEqual(enzymes["NotI"].site, "GCGGCCGC")
        self.assertEqual(enzymes["NotI"].cut_offset, 2)


if __name__ == "__main__":
    unittest.main()
