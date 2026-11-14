"""Site search, sorting, and fragment lengths.

Site search is brute force: a window of the site length slides
along the sequence and is compared with the site. For a circular
plasmid the search runs over the sequence followed by its first
m - 1 bases, so sites spanning the numbering origin are found.

Cut positions are ordered with insertion sort, the sorting
algorithm from class. It runs in O(c^2) for c cuts, but c is
usually under 20 per digest, so the cost is negligible.
"""


def find_sites(seq, site, circular=True):
    """Return the start positions of site in seq (0-based).

    For a circular sequence, sites that wrap around the origin
    are reported with their start near the end of the sequence.
    """
    seq = seq.upper()
    site = site.upper()
    n = len(seq)
    m = len(site)
    positions = []
    if circular:
        # Search seq plus its first m - 1 bases so that a site
        # spanning the origin is found. Only starts in [0, n) count.
        extended = seq + seq[:m - 1]
        for start in range(n):
            if extended[start:start + m] == site:
                positions.append(start)
    else:
        for start in range(n - m + 1):
            if seq[start:start + m] == site:
                positions.append(start)
    return positions


def insertion_sort(values):
    """Return a new list with values sorted (insertion sort)."""
    arr = list(values)
    for j in range(1, len(arr)):
        key = arr[j]
        i = j - 1
        while i >= 0 and arr[i] > key:
            arr[i + 1] = arr[i]
            i = i - 1
        arr[i + 1] = key
    return arr


def cut_positions(plasmid, enzymes):
    """Return the sorted cut positions for a digest, no duplicates.

    enzymes is a list of Enzyme objects. Each recognition site
    gives one cut at (site start + cut offset). For a circular
    plasmid the cut wraps around with mod. Duplicate cuts (two
    enzymes with the same site and offset) are removed so that
    no 0 bp fragment is produced.
    """
    length = plasmid.length
    cuts = []
    for enzyme in enzymes:
        sites = find_sites(plasmid.sequence, enzyme.site, plasmid.circular)
        for start in sites:
            if plasmid.circular:
                cuts.append((start + enzyme.cut_offset) % length)
            else:
                cuts.append(start + enzyme.cut_offset)
    ordered = insertion_sort(cuts)
    deduped = []
    for cut in ordered:
        if len(deduped) == 0 or deduped[-1] != cut:
            deduped.append(cut)
    return deduped


def fragment_lengths(cuts, length, circular=True):
    """Return the fragment lengths (bp) for cut positions.

    - 0 cuts: [] (the plasmid stays circular, reported as uncut).
    - 1 cut: [length] (the circle opens into one linear piece).
    - n cuts: the gaps between consecutive cuts, plus the last
      fragment, which wraps around the origin: (L - c_n) + c_1.
    """
    cuts = insertion_sort(cuts)
    if len(cuts) == 0:
        return []
    if len(cuts) == 1:
        return [length]
    frags = []
    for i in range(len(cuts) - 1):
        frags.append(cuts[i + 1] - cuts[i])
    if circular:
        frags.append((length - cuts[-1]) + cuts[0])
    else:
        frags = [cuts[0]] + frags + [length - cuts[-1]]
    result = []
    for f in frags:
        if f > 0:
            result.append(f)
    return result
