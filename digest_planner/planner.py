"""DigestPlanner: try every digest, rank the ones that distinguish.

For k enzymes there are k single-enzyme digests and k(k-1)/2
two-enzyme digests. Each digest is evaluated on both plasmids,
and the distinguishing ones are ranked by separation score with
a heap (heapq), so the best results come out first.
"""

import heapq

from .compare import distinguishable, separation_score
from .digest import cut_positions, fragment_lengths


class DigestResult:
    """The outcome of one digest on the two plasmids."""

    def __init__(self, scheme, frags_a, frags_b, is_dist, score):
        self.scheme = scheme  # tuple of enzyme names, e.g. ("BamHI", "EcoRI")
        self.frags_a = frags_a
        self.frags_b = frags_b
        self.distinguishable = is_dist
        self.score = score

    def __str__(self):
        return " + ".join(self.scheme) + ": " + str(self.frags_a) \
            + " vs " + str(self.frags_b) + " (score " \
            + str(round(self.score, 2)) + ")"


def format_frags(frags):
    """Format a fragment list for the report."""
    if len(frags) == 0:
        return "uncut"
    parts = []
    for f in sorted(frags):
        parts.append(str(f))
    return ", ".join(parts)


class DigestPlanner:
    """Search all single and double digests for two plasmids."""

    def __init__(self, plasmid_a, plasmid_b, enzymes, threshold=None):
        if plasmid_a.length != plasmid_b.length:
            raise ValueError("Plasmids must have equal length: "
                             + str(plasmid_a.length) + " vs "
                             + str(plasmid_b.length))
        self.plasmid_a = plasmid_a
        self.plasmid_b = plasmid_b
        self.enzymes = enzymes  # dictionary: name -> Enzyme
        self.threshold = threshold  # None for version 1, e.g. 0.10 for v2
        self.results = []

    def enumerate_schemes(self):
        """Every single enzyme and every unordered pair."""
        names = sorted(self.enzymes)
        schemes = []
        for name in names:
            schemes.append((name,))
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                schemes.append((names[i], names[j]))
        return schemes

    def evaluate(self, scheme):
        """Evaluate one digest scheme on both plasmids."""
        scheme_enzymes = []
        for name in scheme:
            scheme_enzymes.append(self.enzymes[name])
        cuts_a = cut_positions(self.plasmid_a, scheme_enzymes)
        cuts_b = cut_positions(self.plasmid_b, scheme_enzymes)
        frags_a = fragment_lengths(cuts_a, self.plasmid_a.length,
                                   self.plasmid_a.circular)
        frags_b = fragment_lengths(cuts_b, self.plasmid_b.length,
                                   self.plasmid_b.circular)
        score = separation_score(frags_a, frags_b)
        is_dist = distinguishable(frags_a, frags_b, self.threshold)
        return DigestResult(scheme, frags_a, frags_b, is_dist, score)

    def run(self):
        """Evaluate every scheme and keep the distinguishing ones."""
        self.results = []
        for scheme in self.enumerate_schemes():
            result = self.evaluate(scheme)
            self.results.append(result)
        return self.results

    def top(self, k=5):
        """Return the k best distinguishing digests, best first.

        heapq is a min-heap, so scores are pushed negated and the
        smallest (most negative, i.e. highest score) pops first.
        """
        heap = []
        for i in range(len(self.results)):
            result = self.results[i]
            if result.distinguishable:
                # The index i breaks ties so that two results with
                # the same score are never compared with each other.
                heapq.heappush(heap, (-result.score, i, result))
        best = []
        while len(heap) > 0 and len(best) < k:
            item = heapq.heappop(heap)
            best.append(item[2])
        return best

    def report(self, k=5):
        """Print the top-k distinguishing digests as a table."""
        best = self.top(k)
        count = 0
        for r in self.results:
            if r.distinguishable:
                count += 1
        if self.threshold is None:
            version = "v1 (exact)"
        else:
            version = "v2 (threshold=" + str(self.threshold) + ")"
        print("Plasmid A: " + str(self.plasmid_a))
        print("Plasmid B: " + str(self.plasmid_b))
        print("Enzymes: " + str(len(self.enzymes))
              + "; digests tried: " + str(len(self.results))
              + "; distinguishing (" + version + "): " + str(count))
        print()
        if len(best) == 0:
            print("No distinguishing digest found.")
            return best
        print("Digest                    Plasmid A (bp)       "
              "Plasmid B (bp)       Score")
        for result in best:
            label = " + ".join(result.scheme)
            print(label.ljust(25) + " "
                  + format_frags(result.frags_a).ljust(20) + " "
                  + format_frags(result.frags_b).ljust(20) + " "
                  + str(round(result.score, 2)))
        return best
