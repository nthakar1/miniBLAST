import math
from Bio.Align import substitution_matrices
BLOSUM62 = substitution_matrices.load("BLOSUM62")

# IMPLEMENTATION NOTES:
"""
UPDATED EXTENSION PIPELINE
Input: 2-hit filtered list of HSSPs
Extension algorithm (previous implementation was more conceptual):
1. Ungapped extension on all HSSPs
2. Filter by S1 threshold
3. Run seed-centered gapped extension on the S1-filtered list
4. Return best of these gapped alignments

REFINED GAPPED EXTENSION
Before, we used a full n x m DP matrix with indexing to access cells. This approach requires us to full the entire matrix, then backtrack from the best cell.
For seed-centered DP, we only compute cells that are visited (instead of computing the entire matrix). Store the 3-levels as dicts using (i, j) as offsets from the seed in (query, ref) since we don't have array-supported global indexing anymore.
For every explored cell (i, j): compute the 3 states M/I/D (if the "parent" cell wasn't computed, store -inf at the current cell instead -- pruning). Assumption is that if the parent wasn't computed, the path to the current cell isn't promising/worth exploring so we cut it off early.
At each layer, check the Xdrop condition (in the first run-through, use Xdrop_gap). A layer represents moving one level outwards from the seed (in a classic DP matrix, we use a source -> sink approach that would start with the seed. Here, we expand in many relative to the seed at once). The explored/computed region ends up looking more like a cloud of cells around the seed rather than a full rectangle.

For backtracking, we can no longer store full backtracking matrices with pointers for each level. Instead, store pointers to the previous (position, state) so that we can reconstruct the aligned strings. Continue traceback until we've gone through all the pointers OR we trigger the Xdrop_gap_final condition Final Xdrop is more aggressive than Xdrop_gap used in extension -- idea is that it is a final check during traceback to make sure we don't have any trailing tails off the alignment that aren't actually good extensions.

PROTEIN ALIGNMENT SUPPORT
For nucleotide alignment, we keep match/mismatch. For protein alignment we want to use standard scores from BLOSUM. To support both, generalized the extension algorithm
"""

def extendFromSeeds(ref, query, seeds, k, s1, matrix, match, mismatch, gapOpen, gapExtend, Xdrop_ungap, Xdrop_gap, Xdrop_final):
    """
    Extends all seeds and returns the best scoring alignment.
    Assuming seeds is given as a tuple of starting indices with (query_kmer_start, ref_kmer_start, seed_score).
    """

    # ungapped on all HSSPs
    ungapped_hits = []
    for q_start, r_start, seed_score in seeds:
        ungapped = ungappedExtension(query, ref, q_start, r_start, seed_score, k, matrix, match, mismatch, Xdrop_ungap)

        # S1-threshold filtering
        if ungapped["score"] >= s1:
            ungapped_hits.append(ungapped)

    #print(f"UNGAPPED HITS: {ungapped_hits}")
    print(f"{len(ungapped_hits)} ungapped seeds passed to gapped extension")

    if not ungapped_hits:
        return None
    
    ungapped_hits.sort(key=lambda x: x["score"], reverse=True)
    total_hits = len(ungapped_hits)
    n_keep = max(1, math.ceil(0.03 * total_hits))
    ungapped_hits = ungapped_hits[:n_keep]
        
    # CHANGED: Sort and keep only the top 3 ungapped hits
    # ungapped_hits.sort(key=lambda x: x["score"], reverse=True)
    # ungapped_hits = ungapped_hits[:3]
    
    best_alignment = None
    best_score = float("-inf")

    # gapped extension for each surviving seed pair
    for hit in ungapped_hits:
        q_seed = hit["q_seed"]
        r_seed = hit["r_seed"]

        gapped = affineGappedExtension(query, ref, q_seed, r_seed, matrix, match, mismatch, gapOpen, gapExtend, Xdrop_gap, Xdrop_final)

        # update best overall alignment
        if gapped["score"] > best_score:
            best_score = gapped["score"]
            best_alignment = gapped

    #print(f"GAPPED: {best_score}")

    return best_alignment

# MAKAYLA UPDATED SEED SCORE TO BE INHERENT TO SEED
def ungappedExtension(query, ref, q_start, r_start, seed_score, k, matrix=None, match=None, mismatch=None, Xdrop=20):
    """
    Ungapped extension with Xdrop centered around the seed.
    Returns dict with score, boundaries, and alignment strings.
    """

    # score the seed 
    # for runtime, now done in bestSeeds, imported as seed_score
    """
    curr_score = 0
    for i in range(k):
        curr_score += scorePair(query[q_start+i], ref[r_start+i], matrix, match, mismatch)
    """

    # extend right
    i = k
    #temp_score = curr_score
    #best_right_score = curr_score
    temp_score = seed_score
    best_right_score = seed_score
    right_end = k

    while q_start + i < len(query) and r_start + i < len(ref):
        temp_score += scorePair(query[q_start+i], ref[r_start+i], matrix, match, mismatch)

        if temp_score > best_right_score:
            best_right_score = temp_score
            right_end = i + 1

        # check Xdrop condition
        if best_right_score - temp_score > Xdrop:
            break

        i += 1
    
    # extend left
    i = 1
    #temp_score = curr_score
    #best_left_score = curr_score
    temp_score = seed_score
    best_left_score = seed_score
    left_start = 0

    while q_start - i >= 0 and r_start - i >= 0:
        temp_score += scorePair(query[q_start-i], ref[r_start-i], matrix, match, mismatch)

        if temp_score > best_left_score:
            best_left_score = temp_score
            left_start = -i

        # check Xdrop condition
        if best_left_score - temp_score > Xdrop:
            break

        i += 1

    q_left = q_start + left_start
    q_right = q_start + right_end
    r_left = r_start + left_start
    r_right = r_start + right_end

    aligned_q = query[q_left:q_right]
    aligned_r = ref[r_left:r_right]

    return {
        "score": best_left_score + best_right_score - seed_score,
        "alignment": (aligned_q, aligned_r),
        "q_range": (q_left, q_right),
        "r_range": (r_left, r_right),
        "q_seed": q_start,
        "r_seed": r_start
    }

# MAKAYLA UPDATED LOOP AT LINE 169
def affineGappedExtension(query, ref, q_seed, r_seed, matrix=None, match=None, mismatch=None, gapOpen=5, gapExtend=2, Xdrop=30, Xdrop_final=100):
    """
    Row-based DP with X-drop extension + affine gap penalties moving rightwards. Returns alignment and associated score. Removed "free ride" edges, since we require the alignment to anchor at the seed.

    States:
    M = match/mismatch layer
    I = insertion in query (gap in reference), gapping row axis
    D = deletion in query (gap in query), gapping col axis
    """

    # 3-level matrix for affine gap support
    M, I, D = {}, {}, {}

    # pointer dict for backtracking
    back = {}

    M[(0, 0)] = scorePair(query[q_seed], ref[r_seed], matrix, match, mismatch)
    best_score = M[(0, 0)]
    best_pos = (0, 0)
    best_state = "M"

    # since we are no longer using the n x m matrix, need to set some limit
    max_band = min(len(query), len(ref)) # limit by length
    # max_band = 2 * (Xdrop // min(match, mismatch)) # limit by scoring scheme

    for d in range(1, max_band):
        layer_best = float("-inf")

        for i in range(-d, d+1):
        #CLAUDE SUGGESTS THIS TO IMPROVE RUNTIME BUT I DONT GET THE LOGIC
        #for i in range(max(-d, -max_band), min(d, max_band) + 1):
            j = d

            qi = q_seed + i
            ri = r_seed + j

            if qi < 0 or qi >= len(query) or ri < 0 or ri >= len(ref):
                continue

            key = (i, j)
            s = scorePair(query[qi], ref[ri], matrix, match, mismatch)
            # print(f"s ({query[qi]}, {ref[ri]}) = {s}")

            prev_diag = (i - 1, j - 1)
            prev_left = (i, j - 1)
            prev_up = (i - 1, j)

            # I state
            i_M = M.get(prev_left, float("-inf")) - gapOpen
            i_I = I.get(prev_left, float("-inf")) - gapExtend

            if i_M >= i_I:
                I[key] = i_M
                back[(key, "I")] = (prev_left, "M")
            else:
                I[key] = i_I
                back[(key, "I")] = (prev_left, "I")

            # D state
            d_M = M.get(prev_up, float("-inf")) - gapOpen
            d_D = D.get(prev_up, float("-inf")) - gapExtend

            if d_M >= d_D:
                D[key] = d_M
                back[(key, "D")] = (prev_up, "M")
            else:
                D[key] = d_D
                back[(key, "D")] = (prev_up, "D")

            # M state
            if prev_diag in M:
                m_M = M[prev_diag] + s
            else:
                m_M = float("-inf")

            candidates = [(m_M, (prev_diag, "M")),
                          (I[key], (key, "I")),
                          (D[key], (key, "D"))]
            best_val, best_back = max(candidates, key=lambda x: x[0])

            M[key] = best_val
            back[(key, "M")] = best_back

            # track best
            cell_score = max(M[key], I[key], D[key])
            # print(f"CELL_SCORE {key}: {cell_score}")

            if cell_score > best_score:
                best_score = cell_score
                best_pos = key

                if cell_score == M[key]:
                    best_state = "M"
                elif cell_score == I[key]:
                    best_state = "I"
                else:
                    best_state = "D"
            
            layer_best = max(layer_best, cell_score)

        # check Xdrop condition
        if best_score - layer_best > Xdrop:
            break

    # print(f"BEST_SCORE: {best_score}")
    # print(f"BEST_POS: {best_pos}")
    # print(f"BEST_STATE: {best_state}")

    # backtracking
    aligned_q, aligned_r = backtrack(query, ref, q_seed, r_seed, back, best_pos, best_state, best_score, Xdrop_final, matrix, match, mismatch, gapOpen, gapExtend)

    q_cov = computeQueryCov(query, aligned_q)
    pct_identity = computePercentIdentity(aligned_q, aligned_r)

    i, j = best_pos
    actual_q_pos = q_seed + i
    actual_r_pos = r_seed + j

    return {
        "score": best_score,
        "alignment": (aligned_q, aligned_r),
        "position": (actual_q_pos, actual_r_pos),
        "query_coverage": q_cov,
        "pct_identity": pct_identity
    }

def backtrack(query, ref, q_seed, r_seed, back, start_pos, start_state, best_score, Xdrop, matrix, match, mismatch, gapOpen, gapExtend):
    """
    Backtracks alignement with affine gap penalties and Xdrop_final condition to check final best alignment and reconstruct aligned strings.
    """

    aligned_q, aligned_r = [], []

    pos = start_pos
    state = start_state

    curr_score = best_score

    while True:
        i, j = pos
        qi = q_seed + i
        ri = r_seed + j

        # add current char to alignment
        # match/mismatch
        if state == "M":
            aligned_q.append(query[qi])
            aligned_r.append(ref[ri])

            score = scorePair(query[qi], ref[ri], matrix, match, mismatch)
        elif state == "I":
            aligned_q.append(query[qi])
            aligned_r.append("-")

            prev_state = back.get((pos, state), (None, None))[1]
            if prev_state == "M":
                score = -gapOpen
            else:
                score = -gapExtend
        elif state == "D":
            aligned_q.append("-")
            aligned_r.append(ref[ri])

            prev_state = back.get((pos, state), (None, None))[1]
            if prev_state == "M":
                score = -gapOpen
            else:
                score = -gapExtend
        else:
            raise ValueError("Unexpected state during backtracking.")
        
        # update running score (working in reverse)
        curr_score += score

        # check Xdrop condition
        if best_score - curr_score > Xdrop:
            break

        # stop if no parent
        if (pos, state) not in back:
            break

        # otherwise update to the prev state
        pos, state = back[(pos, state)]

    aligned_q_str = "".join(reversed(aligned_q))
    aligned_r_str = "".join(reversed(aligned_r))

    return aligned_q_str, aligned_r_str

def computeQueryCov(query, aligned_query):
    """
    Calculates percent query coverage (amount of the input query sequence used in the alignment).
    """

    # find num of chars used in alignment
    aligned_chars = 0
    for symbol in aligned_query:
        if symbol != "-":
            aligned_chars += 1
    
    q_cov = (aligned_chars / len(query)) * 100
    return q_cov

def scorePair(a, b, matrix=None, match=None, mismatch=None):
    """
    General scoring function:
    - protein: takes BLOSUM as substitution matrix
    - nucleotide: no matrix passed, use simple match/mismatch
    """

    if matrix is not None:
        return matrix[a][b]
    else:
        if a == b:
            return match
        else:
            return -mismatch
        
def computeKmerCoverage(query, seeds, k):
    """Returns the percent of query k-mer start positions that produced at least one 
    seed against the reference, calculated as unique seeded positions / total k-mer positions * 100."""
    total_kmer_positions = len(query) - k + 1
    if total_kmer_positions <= 0 or not seeds:
        return 0.0

    seeded_positions = {q_pos for q_pos, _ in seeds}
    valid_seeded = seeded_positions & set(range(total_kmer_positions))
    return (len(valid_seeded) / total_kmer_positions) * 100

def computePercentIdentity(aligned_query, aligned_ref):
    """Returns the percent of alignment positions where both sequences share the 
    same non-gap residue, calculated as identical positions / alignment length * 100."""
    if len(aligned_query) != len(aligned_ref):
        raise ValueError(f"Aligned sequences must be the same length "
                         f"(got {len(aligned_query)} vs {len(aligned_ref)})")
    alignment_length = len(aligned_query)
    if alignment_length == 0:
        return 0.0

    identical = sum(
        1 for a, b in zip(aligned_query, aligned_ref)
        if a == b and a != "-" and b != "-"
    )
    return (identical / alignment_length) * 100

# def main():

#     # query = "ACTGCTG"
#     # ref = "ACTCACTG"
#     query = "MPQLSLSWLGLGPVAASPWLLLLLVGGSWLLARVLAWTYTFYDNCRRLQCFPQPPKQNWFWGHQGLVTPTEEGMKTLTQLVTTYPQGFKLWLGPTFPLLILCHPDIIRPITSASAAVAPKDMIFYGFLKPWLGDGLLLSGGDKWSRHRRMLTPAFHFNILKPYMKIFNKSVNIMHDKWQRLASEGSARLDMFEHISLMTLDSLQKCVFSFESNCQEKPSEYIAAILELSAFVEKRNQQILLHTDFLYYLTPDGQRFRRACHLVHDFTDAVIQERRCTLPTQGIDDFLKNKAKSKTLDFIDVLLLSKDEDGKELSDEDIRAEADTFMFEGHDTTASGLSWVLYHLAKHPEYQEQCRQEVQELLKDREPIEIEWDDLAQLPFLTMCIKESLRLHPPVPVISRCCTQDFVLPDGRVIPKGIVCLINIIGIHYNPTVWPDPEVYDPFRFDQENIKERSPLAFIPFSAGPRNCIGQAFAMAEMKVVLALTLLHFRILPTHTEPRRKPELILRAEGGLWLRVEPLGANSQ"
#     ref = "MSQLSLSWLGLWPVAASPWLLLLLVGASWLLAHVLAWTYAFYDNCRRLRCFPQPPRRNWFWGHQGMVNPTEEGMRVLTQLVATYPQGFKVWMGPISPLLSLCHPDIIRSVINASAAIAPKDKFFYSFLEPWLGDGLLLSAGDKWSRHRRMLTPAFHFNILKPYMKIFNESVNIMHAKWQLLASEGSACLDMFEHISLMTLDSLQKCVFSFDSHCQEKPSEYIAAILELSALVSKRHHEILLHIDFLYYLTPDGQRFRRACRLVHDFTDAVIQERRRTLPSQGVDDFLQAKAKSKTLDFIDVLLLSKDEDGKKLSDEDIRAEADTFMFEGHDTTASGLSWVLYHLAKHPEYQERCRQEVQELLKDREPKEIEWDDLAHLPFLTMCMKESLRLHPPVPVISRHVTQDIVLPDGRVIPKGIICLISVFGTHHNPAVWPDPEVYDPFRFDPENIKERSPLAFIPFSAGPRNCIGQTFAMAEMKVVLALTLLRFRVLPDHTEPRRKPELVLRAEGGLWLRVEPLS"
#     seeds = [(0, 0)]
#     k = 4
#     matrix = BLOSUM62
#     match = None
#     mismatch = None
#     gapOpen = 5
#     gapExtend = 2
#     Xdrop_ungap = 20
#     Xdrop_gap = 30
#     Xdrop_gapFinal = 100
#     s1 = 10

#     #print(BLOSUM62)

#     #print()

#     #print("q:", query)
#     #print("r:", ref)

#     result = extendFromSeeds(ref, query, seeds, k, s1, matrix, match, mismatch, gapOpen, gapExtend, Xdrop_ungap, Xdrop_gap, Xdrop_gapFinal)
    
#     #print(result)



# if __name__ == "__main__":
#     main()