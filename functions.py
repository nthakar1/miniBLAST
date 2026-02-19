# Parent function that returns a maximally scoring local alignment between 2 strings. CURRENTLY INCOMPLETE
# Input: (probably missing things currently)
# Output:
def miniBLAST(ref, query, k, matchScore, mismatchPen, threshHSSP):
    seeds = BestSeeds(ref, query, k, matchScore, mismatchPen, threshHSSP)
    # localAlignments = ExtendFromSeeds(ref, query, seeds, ...)
    pass

# Input: two DNA string query and ref, an int k, scores for matches and mismatches, and a threshold score defining an HSSP
# Output: the highest scoring segment pairs (HSSPs) to be used as seeds for the local alignments. HSSPs will be returned as a dictionary mapping the kmer from the reference to a tuple of indeces where the alignment occured (indexInRef, indexInQuery)
def BestSeeds(ref, query, k, matchScore, mismatchPen, threshold):
    d = Indexation(ref, k)
    HSSPs = {}
    for i in range(len(query)+1-k):
        qKmer = query[i:i+k]
        for rKmer, x in d.items():
            score = ScoreKmers(qKmer, rKmer, matchScore, mismatchPen)
            if score >= threshold:
                HSSPs[rKmer] = (x, i)
    return HSSPs

# Input: a DNA string ref and an int k
# Output: a dictionary mapping each kmer in ref to its start index 
def Indexation(ref, k):
    d = {}
    for i in range(len(ref)+1-k):
        kmer = ref[i:i+k]
        d[kmer] = i
    return d

# Input: two kmers of the same length and two scoring parameters
# Output: the score of the ungapped alignment between the kmers
def ScoreKmers(qKmer, rKmer, matchScore, mismatchPen):
    score = 0
    for i in range(len(qKmer)):
        if qKmer[i] == rKmer[i]:
            score += matchScore
        else:
            score -= mismatchPen
    return score 