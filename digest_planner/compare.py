"""Decide whether two fragment lists would look different on a gel.

Version 1 (baseline): a digest distinguishes the two plasmids
whenever the sorted fragment lists differ at all.

Version 2 (improved): a digest distinguishes them only when at
least one band differs by at least the threshold (a relative
difference, 0.10 by default) from every band in the other lane.
The score is the largest such nearest-band relative difference:

    d(x, Y) = min over y in Y of |x - y| / max(x, y)
    score = the largest d over all bands in both lanes

A score of 0.10 means at least one band differs by 10% or more
from every band in the other lane. Bands a few percent apart
would run as one band on a real gel, so version 2 no longer
reports them as distinguishing.
"""

from .digest import insertion_sort


def differs_exact(a, b):
    """True when the sorted fragment lists differ at all."""
    return insertion_sort(a) != insertion_sort(b)


def nearest_diff(x, others):
    """Smallest relative difference between x and any value in others."""
    best = None
    for y in others:
        diff = abs(x - y) / max(x, y)
        if best is None or diff < best:
            best = diff
    return best


def separation_score(a, b):
    """Largest nearest-band relative difference between two lanes.

    Returns 0.0 when both lanes are identical (including both
    uncut) and 0.0 when exactly one lane is uncut, so that
    one-uncut results rank last in the report.
    """
    if len(a) == 0 or len(b) == 0:
        return 0.0
    score = 0.0
    for x in a:
        score = max(score, nearest_diff(x, b))
    for y in b:
        score = max(score, nearest_diff(y, a))
    return score


def distinguishable(a, b, threshold=None):
    """Decide whether two fragment lists would look different.

    threshold=None (version 1): True when the lists differ.
    threshold set (version 2): True when the separation score
    reaches the threshold. Exactly one plasmid uncut counts as
    distinguishing in both versions, but scores 0.0, so it
    ranks last.
    """
    if (len(a) == 0) != (len(b) == 0):
        return True
    if len(a) == 0 and len(b) == 0:
        return False
    if threshold is None:
        return differs_exact(a, b)
    return separation_score(a, b) >= threshold
