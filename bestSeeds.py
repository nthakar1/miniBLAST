from itertools import product
from extendSeeds import scorePair
from lowComplexityMasker import mask_low_complexity_regions, filter_seeds_by_masking, masking_stats
import time

##### NOTES FOR LATER: #####
"""
- Will want to account for IUPAC nucleotide ambiguity codes during indexation and scoring (reference or query seqs with degenerated bases)
"""
def BestSeedsWithMasking(ref, query, k, matchScore, mismatchPen, matrix, threshHSSP, 
                        apply_masking=True, dust_complexity_threshold=20, mask_char='lowercase', overlap_threshold=0.3):
    """
    BestSeeds with NCBI-standard DUST masking for low complexity regions.
    -----------
    Inputs:
    ref, query (str): DNA or protein strings
    k (int): kmer length
    matchScore, mismatchPen (int): DNA scoring matrix parameters
    matrix (dict or None): BLOSUM-62 substitution matrix
    threshHSSP (float): HSSP threshold for seed scoring
    apply_masking (bool, default True): decide whether to apply masking or not
    dust_complexity_threshold (int, default 20): score threshold that defines complexity of regions. lower scores=higher complexity because less repeats
    mask_char(str): users decide whether they want masked regions to appear as the lowercase version of the original sequence (soft masking) or as a string of Ns (hard masking)
    overlap_threshold(float): determines how much of a seed can overlap with a low complexity region
    -----------    
    Outputs:
    seeds (list): List of seeds with [q_start, r_start, score]. May be filtered by masked sequences or unfiltered depending on selections
    masking_info (dict): Information about masking applied
    """
    # initialize dictionary 
    masking_info = {
        'query_masked': False,
        'ref_masked': False,
        'query_stats': {},
        'ref_stats': {},
        'seeds_before_masking_filter': 0,
        'seeds_after_masking_filter': 0,
        'dust_threshold': dust_complexity_threshold
    }
    
    original_ref = ref
    original_query = query
    query_masked_intervals = []
    ref_masked_intervals = []
    
    # if the apply_masking option was selected
    if apply_masking:
        print("Applying DUST masking for low complexity regions...")
        
        # Mask query sequence (convert low complexity to lowercase), get all the masked intervals
        query, query_masked_intervals = mask_low_complexity_regions(
            query, score_threshold=dust_complexity_threshold, mask_char=mask_char
        )
        masking_info['query_masked'] = len(query_masked_intervals) > 0
        masking_info['query_stats'] = masking_stats(original_query, query_masked_intervals)
        
        # Mask reference sequence  
        ref, ref_masked_intervals = mask_low_complexity_regions(
            ref, score_threshold=dust_complexity_threshold, mask_char=mask_char
        )
        masking_info['ref_masked'] = len(ref_masked_intervals) > 0
        masking_info['ref_stats'] = masking_stats(original_ref, ref_masked_intervals)
        
        if masking_info['query_masked'] or masking_info['ref_masked']:
            print(f"  Query: {masking_info['query_stats']['n_intervals']} intervals, " +
                  f"{masking_info['query_stats']['pct_masked']:.1f}% masked")
            print(f"  Ref:   {masking_info['ref_stats']['n_intervals']} intervals, " +
                  f"{masking_info['ref_stats']['pct_masked']:.1f}% masked")
    
    d = EncodedIndexation(ref, k, matrix)
    seeds = []
    for i in range(len(query)+1-k):
        qKmer = query[i:i+k]

        # generate a list kmers, including qKmer, with an ungapped alignment at or above the HSSP threshold
        potentialHSSPs = GenerateQualifyingHSSPs(qKmer, matchScore, mismatchPen, matrix, threshHSSP)

        
        # if a potential HSSP matches a reference kmer in d, save the qKmer start position (i) and the rKmer start position(s) (d[HSSP])
        for kmer, score in potentialHSSPs:
             check = KmerNumericalEncoding(kmer, matrix)
             if check in d:
                for pos in d[check]:
                    match = [i, pos, score]
                    seeds.append(match)   

    masking_info['seeds_before_masking_filter'] = len(seeds)
    
    # if there are seeds that overlap with low complexity regions AND the user decides to apply masking (default=true), then filter the seeds
    if apply_masking and (query_masked_intervals or ref_masked_intervals):
        seeds = filter_seeds_by_masking(
            seeds, k, query_masked_intervals, ref_masked_intervals, 
            overlap_threshold
        )
    
    masking_info['seeds_after_masking_filter'] = len(seeds)
    
    return seeds, masking_info 
            
# Input: two DNA strings ref and query, an int k, scoring parameters and and HSSP threshold
# Output: a 2d list of seeds represented as [start coord in query, start coord in ref]
def BestSeeds(ref, query, k, matchScore, mismatchPen, matrix, threshHSSP):
    d = EncodedIndexation(ref, k, matrix)
    seeds = []
    for i in range(len(query)+1-k):
        qKmer = query[i:i+k]

        # generate a list kmers, including qKmer, with an ungapped alignment at or above the HSSP threshold
        potentialHSSPs = GenerateQualifyingHSSPs(qKmer, matchScore, mismatchPen, matrix, threshHSSP)

        
        # if a potential HSSP matches a reference kmer in d, save the qKmer start position (i) and the rKmer start position(s) (d[HSSP])
        for kmer, score in potentialHSSPs:
             check = KmerNumericalEncoding(kmer, matrix)
             if check in d:
                for pos in d[check]:
                    match = [i, pos, score]
                    seeds.append(match)  
    return seeds 

    # Best seeds will not contain masking algorithm for now because masking vs no masking logic already included in BestSeedsWithMasking. Using best seeds directly will make protein alignment faster
    # masking_info['seeds_before_masking_filter'] = len(seeds)
    
    # # if there are seeds that overlap with low complexity regions AND the user decides to apply masking (default=true), then filter the seeds
    # if apply_masking and (query_masked_intervals or ref_masked_intervals):
    #     seeds = filter_seeds_by_masking(
    #         seeds, k, query_masked_intervals, ref_masked_intervals, 
    #         overlap_threshold
    #     )
    
    # masking_info['seeds_after_masking_filter'] = len(seeds)
    
    # return seeds, masking_info 

### NEED TO UPDATE THIS TO WORK BETTER FOR IUPAC
# Input: a reference sequence string and an int k
# Output: a dictionary mapping numerically-encoded kmers in ref to a list of their start positions in ref
def EncodedIndexation(ref, k, matrix):
    d = {}
    for i in range(len(ref)+1-k):
        kmer = ref[i:i+k]
        code = KmerNumericalEncoding(kmer, matrix)
        if code is None:
            continue
        if code not in d:
            d[code] = [i]
        else:
            d[code].append(i)
    return d

### NEED TO UPDATE THIS TO WORK BETTER FOR IUPAC
# Input: a kmer string
# Output: a numerical representation of the kmer
def KmerNumericalEncoding(kmer, matrix):
    if matrix == None:
        encodingDict = {"A":0, "C":1, "G":2, "T":3}
    else:
        encodingDict = {"A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7, "K": 8, "L": 9, "M": 10, "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15, "T": 16, "V": 17, "W": 18, "Y": 19}
    
    k = len(kmer)
    res = 0
    for pos in range(k):
        char = kmer[pos]
        # skip if not in encodingDict (ambigious nucleotides, residues)
        if char not in encodingDict:
            return None
        val = encodingDict[char]
        res += val * (4**(k-pos-1))
    return res

# Input: a query kmer and scoring/threshold parameters
# Output: a list of kmers that meet or exceed the HSSP threshold in an ungapped alignment
def GenerateQualifyingHSSPs(qKmer, matchScore, mismatchPen, matrix, threshHSSP):
    # Skip k-mers that are mostly lowercase (masked regions)
    lowercase_count = sum(1 for c in qKmer if c.islower())
    if lowercase_count >= len(qKmer) * 0.7:  # 70% or more masked
        return []  # Return empty list - no candidates
    
    # Convert to uppercase for processing (soft masking approach)
    qKmer = qKmer.upper()

    candidates = []
    scores = []

    while len(candidates) == 0 or len(candidates[0]) < len(qKmer):
        candidates, scores = ExpandCandidates(candidates, scores, qKmer, matrix, matchScore, mismatchPen)
        maxRemainingScore = BestRestOfScore(len(candidates[0]), qKmer, matchScore, matrix)
        keepers = []
        keeper_scores = []
        for x in range(len(candidates)):
            candidate = candidates[x]
            score = scores[x]
            if score + maxRemainingScore >= threshHSSP:
                keepers.append(candidate)
                keeper_scores.append(score)
        candidates = keepers
        scores = keeper_scores

    return [(candidates[i], scores[i]) for i in range(len(candidates))]

def ExpandCandidates(candidates, scores, qKmer, matrix, matchScore, mismatchPen):
    new_candidates = []
    new_scores = []
    alpha_dna = ["A", "C", "T", "G"]
    alpha_aa = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    
    if len(candidates) == 0:
        if matrix == None:
            new_candidates = alpha_dna
        else:
            new_candidates = alpha_aa
        for char in new_candidates:
            score = scorePair(char, qKmer[0], matrix, matchScore, mismatchPen)
            new_scores.append(score)
        return new_candidates, new_scores

    else:
        for i in range(len(candidates)):
            candidate = candidates[i]
            j = len(candidate) # position being added to candidate, to compare to same position in qKmer
            if matrix == None:
                for char in alpha_dna:
                    new_candidate = candidate + char
                    new_score = scores[i] + scorePair(char, qKmer[j], matrix, matchScore, mismatchPen)
                    new_candidates.append(new_candidate)
                    new_scores.append(new_score)
            else:
                for char in alpha_aa:
                    new_candidate = candidate + char
                    new_score = scores[i] + scorePair(char, qKmer[j], matrix, matchScore, mismatchPen)
                    new_candidates.append(new_candidate)
                    new_scores.append(new_score)
    
    return new_candidates, new_scores

def BestRestOfScore(candidateLength, qKmer, matchScore, matrix):
    # for nucleotides, maximum gain is all matches for rest of candidate extension
    if matrix == None:
        lengthLeft = len(qKmer) - candidateLength
        return lengthLeft * matchScore
    
    # for peptides, maximum gains is all matches, score dependent on character
    else:
        seqLeft = qKmer[candidateLength:]
        maxGain = 0
        for char in seqLeft:
            maxGain += matrix[char][char]
        return maxGain

# NEED TO UPDATE THIS TO WORK BETTER FOR IUPAC
# Input: two kmers of the same length and two scoring parameters
# Output: the score of the ungapped alignment between the kmers
def ScoreKmers(qKmer, rKmer, matrix, matchScore, mismatchPen):
    score = 0
    for i in range(len(qKmer)):
        score += scorePair(qKmer[i], rKmer[i], matrix, matchScore, mismatchPen)

    return score 

# Input: an int k
# Output: a list of all possible DNA strings of length k
def GenerateAllKmers(k):
    bases = ["A", "C", "T", "G"]
    return [''.join(p) for p in product(bases, repeat=k)]
def main():
    print("Testing H1N1 sequences with masking...")
    
    start_time = time.perf_counter()
    
    # H1N1 1934 segment 1 (the original long test)
    query = "AGCGAAAGCAGGTCAATTATATTCAATATGGAAAGAATAAAAGAACTAAGAAATCTAATGTCGCAGTCTCGCACCCGCGAGATACTCACAAAAACCACCGTGGACCATATGGCCATAATCAAGAAGTACACATCAGGAAGACAGGAGAAGAACCCAGCACTTAGGATGAAATGGATGATGGCAATGAAATATCCAATTACAGCAGACAAGAGGATAACGGAAATGATTCCTGAGAGAAATGAGCAAGGACAAACTTTATGGAGTAAAATGAATGATGCCGGATCAGACCGAGTGATGGTATCACCTCTGGCTGTGACATGGTGGAATAGGAATGGACCAATGACAAATACAGTTCATTATCCAAAAATCTACAAAACTTATTTTGAAAGAGTCGAAAGGCTAAAGCATGGAACCTTTGGCCCTGTCCATTTTAGAAACCAAGTCAAAATACGTCGGAGAGTTGACATAAATCCTGGTCATGCAGATCTCAGTGCCAAGGAGGCACAGGATGTAATCATGGAAGTTGTTTTCCCTAACGAAGTGGGAGCCAGGATACTAACATCGGAATCGCAACTAACGATAACCAAAGAGAAGAAAGAAGAACTCCAGGATTGCAAAATTTCTCCTTTGATGGTTGCATACATGTTGGAGAGAGAACTGGTCCGCAAAACGAGATTCCTCCCAGTGGCTGGTGGAACAAGCAGTGTGTACATTGAAGTGTTGCATTTGACTCAAGGAACATGCTGGGAACAGATGTATACTCCAGGAGGGGAAGTGAAGAATGATGATGTTGATCAAAGCTTGATTATTGCTGCTAGGAACATAGTGAGAAGAGCTGCAGTATCAGCAGACCCACTAGCATCTTTATTGGAGATGTGCCACAGCACACAGATTGGTGGAATTAGGATGGTAGACATCCTTAAGCAGAACCCAACAGAAGAGCAAGCCGTGGGTATATGCAAGGCTGCAATGGGACTGAGAATTAGCTCATCCTTCAGTTTTGGTGGATTCACATTTAAGAGAACAAGCGGATCATCAGTCAAGAGAGAGGAAGAGGTGCTTACGGGCAATCTTCAAACATTGAAGATAAGAGTGCATGAGGGATATGAAGAGTTCACAATGGTTGGGAGAAGAGCAACAGCCATACTCAGAAAAGCAACCAGGAGATTGATTCAGCTGATAGTGAGTGGGAGAGACGAACAGTCGATTGCCGAAGCAATAATTGTGGCCATGGTATTTTCACAAGAGGATTGTATGATAAAAGCAGTTAGAGGTGATCTGAATTTCGTCAATAGGGCGAATCAGCGACTGAATCCTATGCATCAACTTTTAAGACATTTTCAGAAGGATGCGAAAGTGCTTTTTCAAAATTGGGGAGTTGAACCTATCGACAATGTGATGGGAATGATTGGGATATTGCCCGACATGACTCCAAGCATCGAGATGTCAATGAGAGGAGTGAGAATCAGCAAAATGGGTGTAGATGAGTACTCCAGCACGGAGAGGGTAGTGGTGAGCATTGACCGGTTCTTGAGAGTCCGGGACCAACGAGGAAATGTACTACTGTCTCCCGAGGAGGTCAGTGAAACACAGGGAACAGAGAAACTGACAATAACTTACTCATCGTCAATGATGTGGGAGATTAATGGTCCTGAATCAGTGTTGGTCAATACCTATCAATGGATCATCAGAAACTGGGAAACTGTTAAAATTCAGTGGTCCCAGAACCCTACAATGCTATACAATAAAATGGAATTTGAACCATTTCAGTCTTTAGTACCTAAGGCCATTAGAGGCCAATACAGTGGGTTTGTGAGAACTCTGTTCCAACAAATGAGGGATGTGCTTGGGACATTTGATACCGCACAGATAATAAAACTTCTTCCCTTCGCAGCCGCTCCACCAAAGCAAAGTAGAATGCAGTTCTCCTCATTTACTGTGAATGTGAGGGGATCAGGAATGAGAATACTTGTAAGGGGCAATTCTCCTGTATTCAACTACAACAAGGCCACGAAGAGACTCACAGTTCTCGGAAAGGATGCTGGCACTTTAACCGAAGACCCAGATGAAGGCACAGCTGGAGTGGAGTCCGCTGTTCTGAGGGGATTCCTCATTCTGGGCAAAGAAGACAGGAGATATGGGCCAGCATTAAGCATCAATGAACTGAGCAACCTTGCGAAAGGAGAGAAGGCTAATGTGCTAATTGGGCAAGGAGACGTGGTGTTGGTAATGAAACGAAAACGGGACTCTAGCATACTTACTGACAGCCAGACAGCGACCAAAAGAATTCGGATGGCCATCAATTAGTGTCGAATAGTTTAAAAACGACCTTGTTTCTACT"

    # reference PV175220.1
    ref = "ATGGATGTAAATCCGACTCTACTTTTCTTAAAAGTGCCAGCACAAAATGCTATAAGCACTACATTCCCTTATACTGGAGATCCCCCATATAGCCATGGAACAGGGACAGGATATACCATGGACACAGTCAACAGAACACACCAATACTCAGAAAGGGGCAGATGGACAACAAACACAGAGACGGGAGCACCCCAGCTTAACCCGGTTGATGGACCACTACCTGAGGACAATGAGCCGAGCGGGTATGCTCAAACAGACTGTGTCCTGGAGGCAATGGCTTTTCTTGAAGAATCCCATCCCGGGATATTTGAGAACTCTTGTCTTGAAACGATAGAGGTGGTCCAACAAACAAGAGTGGACAAGCTGACTCAAGGTCGTCAGACCTACGACTGGACATTGAATAGAAACCAACCGGCTGCAACTGCTTTGGCCAATACAATAGAGGTCTTCAGGATGAACAGTCTGACAGCCAATGAGTCGGGAAGGTTAATAGATTTCCTTAAAGATGTAATGGAATCAATGGATAAAGAAGAGATGGAAATAACAACACACTTCCAAAGAAAAAGAAGGATAAGAGACAACATGACCAAGAAAATGGTAACACAAAGAACAATAGGGAAGAAAAAGCATAAACTGAACAAAAGGAGTTATCTAATAAGAGCATTGACACTAAACACGATGACAAAAGATGCAGAAAGAGGCAAATTGAAAAGGCGTGCAATTGCAACACCAGGGATGCAGATCAGAGGATTTGTGTACTTCGTAGAAGCACTGGCGAGAAGCATCTGTGAGAAACTTGAGCAATCTGGACTCCCAGTTGGGGGAAATGAAAAGAAAGCCAAATTGGCAAATGTTGTGAGGAAAATGATGACCAATTCACAAGATACAGAGCTTTCTTTTACAATTACTGGAGACAACACCAAATGGAATGAAAATCAAAACCCCCGAATATTCCTGGCGATGATAACGTACATCACAAGAAATCAGCCTGAATGGTTCAGAAATGTTTTAAGCATTGCCCCTATAATGTTCTCAAACAAAATGGCAAGATTGGGAAAAGGATACATGTTCGAAAGTAAGAGCATGAAGCTGCGAACACAAATACCAGCAGAGATGCTTGCAGATATTGACTTGAAGTACTTCAATGAATCAACAAGGAAAAAGATCGAGAAGATAAGGCCACTCCTCATAGACGGCACAGCCTCATTAAGCCCTGGAATGATGATGGGCATGTTTAACATGCTGAGCACTGTTTTAGGAGTCTCAGTCCTGAACCTCGGGCAAAAGAGATACACTAAAACCACCTACTGGTGGGATGGGCTCCAATCCTCTGATGATTTCGCCCTCATAGTGAATGCACCCGATCATGAAGGAATACAAGCAGGAGTTGATAGGTTTTACAGGACCTGCAAACTGGTTGGAATCAATATGAGCAAAAAGAAATCTTACATAAATAAAACAGGGACATTTGAATTCACAAGCTTTTTCTACCGCTATGGATTTGTAGCTAATTTCAGTATGGAGCTACCCAGCTTTGGGGTATCTGGAGTCAATGAGTCAGCTGACATGAGCATTGGCGTAACAGTAATAAAGAACAACATGATAAACAATGATCTCGGGCCAGCAACAGCCCAAATGGCTCTTCAATTATTCATCAAAGACTACAGGTATACATACCGATGCCACAGAGGTGACATGCAAATCCAAACAAGAAGATCATTTGAGCTAAAAAAACTATGGGAGCAAACCCATTCAAAGGCAGGACTATTGGTTTCGGATGGAGGACCAAATCTGTACAATATCCGGAATCTCCACATCCCAGAAGTCTGCTTGAAATGGGAGCTAATGGATGTAGATTATCGGGGAAGATTGTGTAATCCTCTAAATCCGTTCGTCAACCATAAGGGAATTGAATCCATAAATAGTGCCGTGATAATGCCAACCCATGGTCCGGCTAAAAGCATGGAATATGATGCTGTTGCAACTACACATTCTTGGATTCCAAAAAGGAATCGTTCCATTCTCAATACCAGCCAAAGGGGGGATTCTTGAAGATGAACAAATGTACCAGAAATGCTGCAATCTATTCGAGAAATTCTTCCCTAGCAGTTCATACAGGAGGCCAGTTGGAATTTCAAGCATGGTGGAGGCCATGGTATCTAGGGCCAGAATTGACGCACGGATTGATTTCGAGTCTGGAAGGATGAAGAAAGAAGAATTTACTGAGATCATGAAGATCTGTTCCACCATTGAAGAACTCAGACGGCAAAAGTAG"

    seeds, masking_info = BestSeedsWithMasking(
        ref=ref, query=query, k=7, matchScore=2, mismatchPen=3, matrix=None, threshHSSP=4,
        apply_masking=True, dust_complexity_threshold=20
    )
    
    end_time = time.perf_counter()
    
    print(f"H1N1 Results:")
    print(f"  Seeds found: {len(seeds)}")
    print(f"  Before masking filter: {masking_info['seeds_before_masking_filter']}")
    print(f"  After masking filter: {masking_info['seeds_after_masking_filter']}")
    print(f"  Query masking: {masking_info['query_stats']['pct_masked']:.1f}%")
    print(f"  Ref masking: {masking_info['ref_stats']['pct_masked']:.1f}%")
    print(f"  Time: {end_time - start_time:.4f} seconds")

if __name__ == "__main__":
    main()
