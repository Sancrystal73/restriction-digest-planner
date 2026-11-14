"""Site search, sorting, and fragment-length computation.

Site search is brute force: a window of the site length slides
along the sequence and is compared with the site. For a circular
plasmid the search runs over the sequence followed by its first
m - 1 bases, so sites spanning the numbering origin are found.

Cut positions are ordered with selection sort (O(c^2) for c cuts;
c is usually under 20 per digest, so the cost is negligible).
"""


def find_sites(seq, site, circular=True):
    """Return 0-based start positions of *site* in *seq*.

    For a circular sequence, sites that wrap around the origin
    are reported with their start position near the end of the
    sequence. Overlapping sites are all reported.
    """
    seq = seq.upper()
    site = site.upper()
    n = len(seq)
    m = len(site)
    if m == 0 or n == 0 or m > n:
        # A site longer than the plasmid cannot occur, even wrapping,
        # unless the plasmid itself is shorter than the site.
        if m > n and circular:
            # Still impossible: a site needs m distinct positions.
            return []
        if m == 0 or n == 0:
            return []
    positions = []
    if circular:
        # Search seq + first (m - 1) bases so wrap-around sites are found.
        # Only starts in [0, n) are valid.
        extended = seq + seq[: m - 1]
        for start in range(n):
            match = True
            for j in range(m):
                if extended[start + j] != site[j]:
                    match = False
                    break
            if match:
                positions.append(start)
    else:
        for start in range(n - m + 1):
            match = True
            for j in range(m):
                if seq[start + j] != site[j]:
                    match = False
                    break
            if match:
                positions.append(start)
    return positions


def selection_sort(values):
    """Return a new list with *values* sorted (selection sort)."""
    arr = list(values)
    n = len(arr)
    for i in range(n):
        min_idx = i
        for j in range(i + 1, n):
            if arr[j] < arr[min_idx]:
                min_idx = j
        if min_idx != i:
            arr[i], arr[min_idx] = arr[min_idx], arr[i]
    return arr


def cut_positions(plasmid, enzymes):
    """Return sorted, de-duplicated cut positions for a digest.

    *enzymes* is an iterable of Enzyme objects. Each recognition
    site contributes one cut at (site_start + cut_offset) mod L
    for a circular plasmid. Duplicate cut positions (e.g. two
    enzymes sharing a site and offset) are removed so no 0 bp
    fragment is produced.
    """
    length = plasmid.length
    cuts = []
    for enzyme in enzymes:
        for start in find_sites(plasmid.sequence, enzyme.site, circular=plasmid.circular):
            if plasmid.circular:
                cut = (start + enzyme.cut_offset) % length
            else:
                cut = start + enzyme.cut_offset
            cuts.append(cut)
    ordered = selection_sort(cuts)
    # Remove duplicates, preserving order.
    deduped = []
    for cut in ordered:
        if not deduped or deduped[-1] != cut:
            deduped.append(cut)
    return deduped


def fragment_lengths(cuts, length, circular=True):
    """Return fragment lengths (bp) for sorted cut positions.

    - 0 cuts: [] (the plasmid stays circular, reported as uncut).
    - 1 cut: [length] (the circle is opened into one linear piece).
    - n cuts: differences between consecutive cuts, plus the
      wrap-around fragment (L - c_n) + c_1 for circular DNA.

    Fragments always sum to *length* when there is at least one cut.
    """
    cuts = selection_sort(cuts)
    if not cuts:
        return []
    if len(cuts) == 1:
        return [length]
    frags = []
    for i in range(len(cuts) - 1):
        frags.append(cuts[i + 1] - cuts[i])
    if circular:
        frags.append((length - cuts[-1]) + cuts[0])
    else:
        # Linear DNA: end pieces are fragments too.
        frags = [cuts[0]] + frags + [length - cuts[-1]]
    # Drop any zero-length fragment (should not happen after dedup).
    return [f for f in frags if f > 0]
