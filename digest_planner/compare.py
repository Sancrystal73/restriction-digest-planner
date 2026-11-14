"""Distinguishability tests and scoring for fragment lists.

Version 1 (baseline): a digest distinguishes the two plasmids
whenever the sorted fragment lists differ at all.

Version 2 (improved): a digest distinguishes them only when at
least one band differs by at least *threshold* (a relative
difference, 0.10 by default) from every band in the other lane.
The score is the largest such nearest-band relative difference:

    d(x, Y) = min_y |x - y| / max(x, y)
    score = max(max_a d(a, B), max_b d(b, A))

A score of 0.10 means at least one band differs by 10% or more
from every band in the other lane. Fragments differing by only a
few percent would run as one band on a real gel, so version 2
no longer reports them as distinguishing.
"""

from .digest import selection_sort


def differs_exact(a, b):
    """True when the sorted fragment lists differ at all."""
    return selection_sort(a) != selection_sort(b)


def _nearest_relative_diff(x, others):
    """Smallest relative difference between x and any value in others."""
    best = None
    for y in others:
        denom = max(x, y)
        diff = abs(x - y) / denom if denom else 0.0
        if best is None or diff < best:
            best = diff
    return best if best is not None else 0.0


def separation_score(a, b):
    """Largest nearest-band relative difference between two lanes.

    Returns 0.0 when both lanes are identical (including both
    uncut) and 0.0 when exactly one lane is uncut, so that
    one-uncut results rank last in the report.
    """
    if not a and not b:
        return 0.0
    if not a or not b:
        return 0.0
    score = 0.0
    for x in a:
        score = max(score, _nearest_relative_diff(x, b))
    for y in b:
        score = max(score, _nearest_relative_diff(y, a))
    return score


def distinguishable(a, b, threshold=None):
    """Decide whether two fragment lists would look different.

    - threshold=None (version 1): True when the sorted lists differ.
    - threshold set (version 2): True when separation_score(a, b)
      is at least the threshold.
    - Exactly one plasmid uncut: True in both versions (the lanes
      differ), but the score is 0.0 so it ranks last.
    """
    if (not a) != (not b):
        return True
    if not a and not b:
        return False
    if threshold is None:
        return differs_exact(a, b)
    return separation_score(a, b) >= threshold
