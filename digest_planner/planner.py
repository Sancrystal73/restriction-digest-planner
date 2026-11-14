"""DigestPlanner: enumerate, evaluate, rank, and report digests.

For k enzymes there are k single-enzyme digests and k(k-1)/2
two-enzyme digests. Each digest is evaluated on both plasmids,
and distinguishing digests are ranked by separation score with
a heap (heapq), so the k best are retrieved without fully
sorting every result.
"""

import heapq
from dataclasses import dataclass, field

from .compare import distinguishable, separation_score
from .digest import cut_positions, fragment_lengths


@dataclass(order=False)
class DigestResult:
    """Outcome of one digest scheme on the two plasmids."""

    scheme: tuple  # enzyme names, e.g. ("EcoRI", "BamHI")
    frags_a: list = field(default_factory=list)
    frags_b: list = field(default_factory=list)
    distinguishable: bool = False
    score: float = 0.0

    @property
    def label(self):
        return " + ".join(self.scheme)

    @property
    def one_uncut(self):
        return (not self.frags_a) != (not self.frags_b)


def _fmt_frags(frags):
    if not frags:
        return "uncut"
    return ", ".join(str(f) for f in sorted(frags))


class DigestPlanner:
    """Search all single and double digests for two plasmids."""

    def __init__(self, plasmid_a, plasmid_b, enzymes, threshold=None):
        """
        Args:
            plasmid_a, plasmid_b: Plasmid objects of equal length.
            enzymes: dict of name -> Enzyme.
            threshold: None for version 1 (exact difference), or a
                float such as 0.10 for version 2.
        """
        if plasmid_a.length != plasmid_b.length:
            raise ValueError(
                f"Plasmids must have equal length: {plasmid_a.length} vs {plasmid_b.length}"
            )
        self.plasmid_a = plasmid_a
        self.plasmid_b = plasmid_b
        self.enzymes = dict(enzymes)
        self.threshold = threshold
        # Cache: (plasmid id, enzyme name) -> cut positions.
        # Lets every pair digest reuse the single-enzyme searches.
        self._cut_cache = {}
        self.results = []
        self._heap = []

    def _cuts_for(self, plasmid, enzyme):
        key = (id(plasmid), enzyme.name)
        if key not in self._cut_cache:
            self._cut_cache[key] = cut_positions(plasmid, [enzyme])
        return self._cut_cache[key]

    def _frags_for(self, plasmid, scheme_enzymes):
        # Merge cached single-enzyme cuts, sort, dedup, then fragment.
        from .digest import selection_sort

        merged = []
        for enzyme in scheme_enzymes:
            merged.extend(self._cuts_for(plasmid, enzyme))
        ordered = selection_sort(merged)
        deduped = []
        for cut in ordered:
            if not deduped or deduped[-1] != cut:
                deduped.append(cut)
        return fragment_lengths(deduped, plasmid.length, circular=plasmid.circular)

    def enumerate_schemes(self):
        """Every single enzyme and every unordered pair."""
        names = sorted(self.enzymes)
        schemes = [(name,) for name in names]
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                schemes.append((names[i], names[j]))
        return schemes

    def evaluate(self, scheme):
        """Evaluate one digest scheme on both plasmids."""
        scheme_enzymes = [self.enzymes[name] for name in scheme]
        frags_a = self._frags_for(self.plasmid_a, scheme_enzymes)
        frags_b = self._frags_for(self.plasmid_b, scheme_enzymes)
        score = separation_score(frags_a, frags_b)
        is_dist = distinguishable(frags_a, frags_b, threshold=self.threshold)
        return DigestResult(
            scheme=tuple(scheme),
            frags_a=frags_a,
            frags_b=frags_b,
            distinguishable=is_dist,
            score=score,
        )

    def run(self):
        """Evaluate every scheme; keep distinguishing ones on a heap.

        Returns the list of all DigestResults. The heap holds
        (-score, index, result) so the highest score pops first.
        """
        self.results = []
        self._heap = []
        for idx, scheme in enumerate(self.enumerate_schemes()):
            result = self.evaluate(scheme)
            self.results.append(result)
            if result.distinguishable:
                heapq.heappush(self._heap, (-result.score, idx, result))
        return self.results

    def top(self, k=5):
        """Return the k best distinguishing digests, best first."""
        if not self._heap and not self.results:
            self.run()
        # Pop from a copy so repeated calls still work.
        heap_copy = list(self._heap)
        heapq.heapify(heap_copy)
        best = []
        while heap_copy and len(best) < k:
            _, _, result = heapq.heappop(heap_copy)
            best.append(result)
        return best

    def report(self, k=5):
        """Print the top-k distinguishing digests as a plain-text table."""
        best = self.top(k)
        total = len(self.results) if self.results else len(self.enumerate_schemes())
        n_dist = sum(1 for r in self.results if r.distinguishable) if self.results else len(best)
        version = "v1 (exact)" if self.threshold is None else f"v2 (threshold={self.threshold})"
        print(f"Plasmid A: {self.plasmid_a.name} ({self.plasmid_a.length} bp)")
        print(f"Plasmid B: {self.plasmid_b.name} ({self.plasmid_b.length} bp)")
        print(f"Enzymes: {len(self.enzymes)}; digests tried: {total}; "
              f"distinguishing ({version}): {n_dist}")
        print()
        if not best:
            print("No distinguishing digest found.")
            return best
        print(f"{'Digest':<25} {'Plasmid A (bp)':<20} {'Plasmid B (bp)':<20} Score")
        for result in best:
            print(f"{result.label:<25} {_fmt_frags(result.frags_a):<20} "
                  f"{_fmt_frags(result.frags_b):<20} {result.score:.2f}")
        one_uncut = [r for r in best if r.one_uncut]
        if one_uncut:
            print()
            print("Note: digests where exactly one plasmid is uncut are listed")
            print("separately in effect (score 0.00, ranked last) because uncut")
            print("plasmid does not run at its true size on a gel.")
        return best
