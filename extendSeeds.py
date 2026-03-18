# IMPLEMENTATION NOTES:
"""
Implementing affine gapped X-drop extension algorithm here. Similar to Smith-Waterman (local alignment recurrence from class), but focusing on a smaller section of the DP matrix.

Assumes correct generation of HSSPs with starting indices -- for extension, we extend outward bidirectionally (anchored in the center at the HSSP).
Uses DP recurrence with affine gap penalties (3-layer DP matrix) to allow indels. X-drop heuristic pruns paths that score more than X below the best score observed so far (allows us to focus on
more promising alignment regions). So keep track of the best global score, then drop threshold becomes (best_global_score - X).

Affine gaps from class vs. in this X-drop version:
Function from class uses a 3-layer DP full matrix approach to do global alignment with affine gaps -- compute every cell, O(nm) space. X-drop version will focus only on the seeded region, so we
use row-based instead -- keep track of curr and prev row since recurrence only depends on direct neighbors, O(m).
Using Smith-Waterman-like recurrence for local alignment (alignment can start anywhere in the sequences). Since class function is global alignment, it uses best score = sink node. For local 
alignment, best score can end anywhere and needs to be updated to keep track of the global best -- this is where we start backtracking from.
Keeping the same affine gap logic for moving between layers + X-drop pruning rule: after computing cell score, check that the score hasn't dropped too low (if cell score > best score - Xdrop,
continue. otherwise prune/stop extending). Terminate the row early if there are no good cells left (instead of completing the full row/matrix anyways in traditional DP approach).
Key point: original asks for the best alignment between two sequences. X-drop asks: if we already have a good matching region/seed, how far can we extend before alignment becomes bad?

Cell-based vs Row-based X-drop:
Found references to a row-based heuristic for X-drop pruning, where the idea is to keep track of the best score in a row and to prune if best_score - row_best > xdrop. we assume that if even
the best cell in the row is too far below the best alignment so far, then we will not be able to recover a good alignment later in the row and should terminate now. Cell-based is more expensive,
so row-based pruning heuristic prunes an entire frontier at once.
Real BLAST uses cell-based. Allows for more irregularly shaped alignments -- useful in cases when most cells have bad scores but one diagonal had recoverable alignment later on (large indels,
temp low-similarity regions/mutation hotspots?, frameshift-like gaps? -- probably more common AA alignment).
Ex: query: AAAAAAAA--------CCCCCCCC
    ref:   AAAAAAAAGGGGGGGGCCCCCCCC
          [--HSSP--]      [--HSSP--]
Depending on chosen params, row-based would drop off early due to G chain indel but cell-based can recover the alignment at the C chain.

Params needed (+ BLAST defaults):
- X (20), can test 10-30 for nucleotides, 40-100 for AAs?
- Match reward
- Mismatch penalty
- Gap opening penalty
- Gap extension penalty

Using row-by-row DP matrix for now (like what we did in class). Pruning with X-drop can shrink the search space quickly, but we'll still be iterating across cols for each row (even if we do
nothing with that cell).
Optimization idea: BLAST uses diagonal bands (diag = q_idx - r_idx) -- idea is that good alignment mostly stay on the diagonal (same as in 2-hit seeding), even with gaps we shouldn't move too
much from the line. Restrict computation in DP matrix to only cells near the seed diagonal (band width depends on X param).
"""

def extendFromSeeds(ref, query, seeds, k, matchReward, mismatchPenalty, gapOpenPenalty, gapExtendPenalty, Xdrop):
    """
    Extends all seeds and returns the best scoring alignment.
    Assuming seeds is given as a tuples of starting indices with (query_kmer_start, ref_kmer_start).
    """

    best_alignment = None
    best_score = float("-inf")

    for q_start, r_start in seeds:
        q_seed = query[q_start:q_start+k]
        r_seed = ref[r_start:r_start+k]

        seed_score = scoreKmers(q_seed, r_seed, matchReward, mismatchPenalty)

        # extend rightwards
        q_right = query[q_start+k:]
        r_right = ref[r_start+k:]

        right = affineExtension_perCell(q_right, r_right, matchReward, mismatchPenalty, gapOpenPenalty, gapExtendPenalty, Xdrop)

        # extend leftwards (using reversed seq)
        q_left = query[:q_start][::-1]
        r_left = ref[:r_start][::-1]

        left = affineExtension_perCell(q_left, r_left, matchReward, mismatchPenalty, gapOpenPenalty, gapExtendPenalty, Xdrop)

        # combine extensions into a single scored alignment
        score = left["score"] + seed_score + right["score"]
        alignment = combineAlignments(left, q_seed, r_seed, right)

        # update best scoring alignment
        if score > best_score:
            best_score = score
            best_alignment = {
                "score": score,
                "query": alignment[0],
                "ref": alignment[1]
            }

    return best_alignment


def combineAlignments(left, q_seed, r_seed, right):
    """
    Combines left extension + seed + right extension.
    """
    q_left, r_left = left["alignment"]
    q_right, r_right = right["alignment"]

    # restoring correct order for left extension
    q_left, r_left = q_left[::-1], r_left[::-1]

    q_aligned = q_left + q_seed + q_right
    r_aligned = r_left + r_seed + r_right
    
    return q_aligned, r_aligned

def affineExtension_perRow(query, ref,
                           matchReward, mismatchPenalty,
                           gapOpenPenalty, gapExtendPenalty,
                           Xdrop):
    """
    Row-based DP with X-drop extension + affine gap penalties moving rightwards. Returns alignment and associated score.

    States:
    M = match/mismatch layer
    I = insertion in query (gap in reference), gapping row axis
    D = deletion in query (gap in query), gapping col axis
    """

    n = len(ref)
    m = len(query)

    NEG_INF = float("-inf")

    # DP matrices
    M = [[0]*(m+1) for _ in range(n+1)]
    I = [[NEG_INF]*(m+1) for _ in range(n+1)]
    D = [[NEG_INF]*(m+1) for _ in range(n+1)]

    # backtracking pointers: store (prev_state, prev_i, prev_j)
    ptr_M = [[None]*(m+1) for _ in range(n+1)]
    ptr_I = [[None]*(m+1) for _ in range(n+1)]
    ptr_D = [[None]*(m+1) for _ in range(n+1)]

    best_score = 0
    best_pos = (0, 0)
    best_state = "M"

    for i in range(1, n+1):

        row_best = NEG_INF

        for j in range(1, m+1):

            # match / mismatch
            if ref[i-1] == query[j-1]:
                diag = matchReward
            else:
                diag = -mismatchPenalty

            # i state -- gap in ref
            from_M = M[i][j-1] - gapOpenPenalty
            from_I = I[i][j-1] - gapExtendPenalty

            if from_M >= from_I:
                I[i][j] = from_M
                ptr_I[i][j] = ("M", i, j-1)
            else:
                I[i][j] = from_I
                ptr_I[i][j] = ("I", i, j-1)

            # d state -- gap in query
            from_M = M[i-1][j] - gapOpenPenalty
            from_D = D[i-1][j] - gapExtendPenalty

            if from_M >= from_D:
                D[i][j] = from_M
                ptr_D[i][j] = ("M", i-1, j)
            else:
                D[i][j] = from_D
                ptr_D[i][j] = ("D", i-1, j)

            # m state -- match/mismatch, pull from both query and ref
            candidates = [
                (0, None),  # restart
                (M[i-1][j-1] + diag, ("M", i-1, j-1)),
                (I[i][j], ("I", i, j)),
                (D[i][j], ("D", i, j))
            ]

            best_val, best_ptr = max(candidates, key=lambda x: x[0])

            M[i][j] = best_val
            ptr_M[i][j] = best_ptr

            # track best overall
            cell_score = max(M[i][j], I[i][j], D[i][j])

            if cell_score > best_score:
                best_score = cell_score
                best_pos = (i, j)

                # track which state we ended in
                if cell_score == M[i][j]:
                    best_state = "M"
                elif cell_score == I[i][j]:
                    best_state = "I"
                else:
                    best_state = "D"

            if cell_score > row_best:
                row_best = cell_score

        # check X-drop pruning condition
        if best_score - row_best > Xdrop:
            break

    # backtracking
    aligned_q = []
    aligned_r = []

    i, j = best_pos
    state = best_state

    while i > 0 and j > 0:

        if state == "M":
            ptr = ptr_M[i][j]
            if ptr is None:
                break

            prev_state, pi, pj = ptr

            if M[i][j] == 0:
                break

            aligned_q.append(query[j-1])
            aligned_r.append(ref[i-1])

        elif state == "I":
            prev_state, pi, pj = ptr_I[i][j]
            aligned_q.append(query[j-1])
            aligned_r.append("-")

        elif state == "D":
            prev_state, pi, pj = ptr_D[i][j]
            aligned_q.append("-")
            aligned_r.append(ref[i-1])

        i, j = pi, pj
        state = prev_state

    aligned_q = "".join(reversed(aligned_q))
    aligned_r = "".join(reversed(aligned_r))

    return {
        "score": best_score,
        "alignment": (aligned_q, aligned_r)
    }

def affineExtension_perCell(query, ref, matchReward, mismatchPenalty, gapOpenPenalty, gapExtendPenalty, Xdrop):
    """
    Row-based DP with X-drop extension + affine gap penalties moving rightwards. Returns alignment and associated score.

    States:
    M = match/mismatch layer
    I = insertion in query (gap in reference), gapping row axis
    D = deletion in query (gap in query), gapping col axis
    """
    n = len(ref)
    m = len(query)
    
    # prev row (don't need I_prev since we move horizontally in the same row)
    M_prev = [0] * (m + 1)
    D_prev = [float("-inf")] * (m + 1)

    best_score = 0
    best_position = (0, 0)

    # iterate through rows (ref string)
    for i in range(1, n + 1):
        # curr row
        M_curr = [0] * (m + 1)
        I_curr = [float("-inf")] * (m + 1)
        D_curr = [float("-inf")] * (m + 1)

        active_cells = False # for tracking whether any cell survives pruning

        # iterate through cols (query string)
        for j in range(1, m + 1):
            # match/mismatch score
            if ref[i-1] == query[j-1]:
                diag = matchReward
            else:
                diag = -mismatchPenalty
            
            # insertion (gap in ref)
            I_curr[j] = max(
                M_curr[j-1] - gapOpenPenalty,
                I_curr[j-1] - gapExtendPenalty
            )

            # deletion (gap in query)
            D_curr[j] = max(
                M_prev[j] - gapOpenPenalty,
                D_prev[j] - gapExtendPenalty
            )

            # find best score for curr cell
            M_curr[j] = max(
                0,
                M_prev[j-1] + diag,
                I_curr[j],
                D_curr[j]
            )

            # best state at this cell
            cell_score = max(M_curr[j], I_curr[j], D_curr[j])
            
            # Xdrop pruning
            if best_score - cell_score > Xdrop:
                M_curr[j] = float("-inf")
                I_curr[j] = float("-inf")
                D_curr[j] = float("-inf")
                continue

            active_cells = True

            # find best global score
            if M_curr[j] > best_score:
                best_score = M_curr[j]
                best_position = (i,j)

        # if entire row is pruned, stop extending
        if not active_cells:
            break

        # slide windows
        M_prev = M_curr
        D_prev = D_curr

    return best_score, best_position

def scoreKmers(qKmer, rKmer, matchReward, mismatchPenalty):
    """
    Scores ungapped alignment between two kmers (LCS score).
    Using simple equality for nucleotides for now, will have to change logic for AA sequences later based on BLOSUM.
    """

    score = 0
    for q, r in zip(qKmer, rKmer):
        if q == r:
            score += matchReward
        else:
            score -= mismatchPenalty
    return score