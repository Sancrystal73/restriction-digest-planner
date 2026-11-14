"""End-to-end tests for planner.py: T13, T14, T16 (optional)."""

import unittest
from pathlib import Path

from digest_planner.enzymes import load_enzymes
from digest_planner.planner import DigestPlanner
from digest_planner.sequence_io import Plasmid

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def load_synthetic():
    enzymes = load_enzymes(DATA / "enzymes.csv")
    a = Plasmid.from_fasta(DATA / "plasmid_A.fasta")
    b = Plasmid.from_fasta(DATA / "plasmid_B.fasta")
    return a, b, enzymes


class TestEnumeration(unittest.TestCase):
    def test_scheme_count(self):
        a, b, enzymes = load_synthetic()
        planner = DigestPlanner(a, b, enzymes)
        schemes = planner.enumerate_schemes()
        k = len(enzymes)
        self.assertEqual(len(schemes), k + k * (k - 1) // 2)  # 12 + 66 = 78
        self.assertIn(("EcoRI",), schemes)
        self.assertIn(("BamHI", "EcoRI"), schemes)


class TestSyntheticPair(unittest.TestCase):
    def test_t13_end_to_end(self):
        a, b, enzymes = load_synthetic()
        planner = DigestPlanner(a, b, enzymes, threshold=0.10)
        planner.run()

        # Single digests are not distinguishing.
        single_ecori = planner.evaluate(("EcoRI",))
        self.assertEqual(single_ecori.frags_a, [3000])
        self.assertEqual(single_ecori.frags_b, [3000])
        self.assertFalse(single_ecori.distinguishable)

        single_bamhi = planner.evaluate(("BamHI",))
        self.assertEqual(single_bamhi.frags_a, [3000])
        self.assertEqual(single_bamhi.frags_b, [3000])
        self.assertFalse(single_bamhi.distinguishable)

        # Double digest distinguishes, with the worked-example fragments.
        double = planner.evaluate(("BamHI", "EcoRI"))
        self.assertEqual(sorted(double.frags_a), [800, 2200])
        self.assertEqual(sorted(double.frags_b), [1400, 1600])
        self.assertTrue(double.distinguishable)
        self.assertAlmostEqual(double.score, 0.43, places=2)

        # It is the only distinguishing digest, and top(1) returns it.
        best = planner.top(1)
        self.assertEqual(len(best), 1)
        self.assertEqual(best[0].scheme, ("BamHI", "EcoRI"))

    def test_top_respects_k(self):
        a, b, enzymes = load_synthetic()
        planner = DigestPlanner(a, b, enzymes, threshold=None)
        planner.run()
        self.assertLessEqual(len(planner.top(5)), 5)


class TestOrientationPair(unittest.TestCase):
    @unittest.skipUnless(
        (DATA / "puc19_forward.fasta").exists()
        and (DATA / "puc19_reverse.fasta").exists(),
        "Orientation pair not generated: run data/make_orientation_pair.py "
        "after downloading pUC19 (L09137) and lambda (J02459).",
    )
    def test_t14_end_to_end(self):
        enzymes = load_enzymes(DATA / "enzymes.csv")
        a = Plasmid.from_fasta(DATA / "puc19_forward.fasta")
        b = Plasmid.from_fasta(DATA / "puc19_reverse.fasta")
        self.assertEqual(a.length, 3686)
        self.assertEqual(b.length, 3686)
        planner = DigestPlanner(a, b, enzymes, threshold=0.10)
        planner.run()

        # BamHI alone cuts both into 1,000 + 2,686 bp: not distinguishing.
        bamhi = planner.evaluate(("BamHI",))
        self.assertEqual(sorted(bamhi.frags_a), [1000, 2686])
        self.assertEqual(sorted(bamhi.frags_b), [1000, 2686])
        self.assertFalse(bamhi.distinguishable)

        # At least one double digest should distinguish the orientations.
        self.assertTrue(any(r.distinguishable for r in planner.results))


class TestBiopythonCrossCheck(unittest.TestCase):
    @unittest.skipUnless(
        (DATA / "puc19.fasta").exists(), "pUC19 FASTA not downloaded."
    )
    def test_t16_optional_biopython(self):
        try:
            from Bio import Restriction
            from Bio.Seq import Seq
        except ImportError:
            self.skipTest("Biopython not installed (optional check).")
        from digest_planner.digest import cut_positions

        enzymes = load_enzymes(DATA / "enzymes.csv")
        puc19 = Plasmid.from_fasta(DATA / "puc19.fasta")
        for name, enzyme in enzymes.items():
            batch = Restriction.RestrictionBatch([name])
            # Bio.Restriction cut positions are 1-based; convert.
            expected = sorted(
                (pos - 1) % puc19.length
                for pos in batch.search(Seq(puc19.sequence), linear=False)[0]
            )
            got = cut_positions(puc19, [enzyme])
            self.assertEqual(got, expected, f"Mismatch for {name}")


if __name__ == "__main__":
    unittest.main()
