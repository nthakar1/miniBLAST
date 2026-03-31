def TwoHitSeeds(single_hits, k, A):
    """
    Implement the 2-Hit Seeding Approach for sequence alignment.
    
    Inputs:
    - ref, query: DNA strings
    - k: word/kmer length
    - matchScore, mismatchPen: scoring parameters
    - T: Score threshold for a single word hit (should be lower than 1-hit threshold)
    - A: Maximum window distance between two hits on the same diagonal
    
    Outputs:
    - A 2D list of valid seeds [start coord in query, start coord in ref] that triggered the 2-hit rule.
    """
    
    # scan for all single hits using the lower threshold T\
    """
    d = EncodedIndexation(ref, k)
    single_hits = []
    
    for i in range(len(query) + 1 - k):
        qKmer = query[i:i+k]
        
        # generate a list of kmers with an ungapped alignment at or above threshold T
        allPossibleKmers = GenerateAllKmers(k)
        potentialHSSPs = []
        for kmer in allPossibleKmers:
            score = ScoreKmers(qKmer, kmer, matchScore, mismatchPen)
            if score >= T:
                potentialHSSPs.append(kmer)

        # check against reference dictionary
        for kmer in potentialHSSPs:
            check = KmerNumericalEncoding(kmer)
            if check in d:
                for pos in d[check]:
                    single_hits.append([i, pos]) # [query_start, ref_start]
    """

    # group hits by their diagonal (q_pos - r_pos)
    diagonals = {}
    for q_pos, r_pos, score in single_hits:
        diag = q_pos - r_pos
        if diag not in diagonals:
            diagonals[diag] = []
        diagonals[diag].append([q_pos, r_pos, score])

    # check proximity and trigger extension
    valid_seeds = []
    
    for diag, hits in diagonals.items():
        # Sort sequentially by query coordinate to scan left-to-right
        hits.sort(key=lambda x: x[0])
        
        last_hit = None
        for curr_hit in hits:
            # Initialize the first hit on this diagonal
            if last_hit is None:
                last_hit = curr_hit
                continue
            
            # Distance between their first coordinates (query positions)
            distance = curr_hit[0] - last_hit[0]
            
            # Condition: Non-overlapping AND within window distance A
            if distance >= k and distance <= A:
                # Valid 2-hit found! 
                # We append the current hit as the seed to trigger bidirectional extension
                valid_seeds.append(curr_hit)
                
                # Update last_hit to curr_hit so we can chain subsequent hits if they occur
                last_hit = curr_hit 
            
            elif distance > A:
                # Hits are too far apart. current hit resets as the new baseline
                last_hit = curr_hit
                
            # If distance < k, they overlap. just ignore curr_hit and keep last_hit.

    return valid_seeds
