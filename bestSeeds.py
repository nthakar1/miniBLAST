from itertools import product
from extendSeeds import scorePair

##### NOTES FOR LATER: #####
"""
- Will want to account for IUPAC nucleotide ambiguity codes during indexation and scoring (reference or query seqs with degenerated bases)
"""

"""
# Input: two DNA string query and ref, an int k, scores for matches and mismatches, and a threshold score defining an HSSP
# Output: the highest scoring segment pairs (HSSPs) to be used as seeds for the local alignments. 
def BestSeeds(ref, query, k, matchScore, mismatchPen, threshHSSP):
    d = Indexation(ref, k)
    HSSPs = {} # a dictionary mapping the kmer from the reference to a tuple of indeces where the alignment occured ([indecesInRef], indexInQuery)
    for i in range(len(query)+1-k):
        qKmer = query[i:i+k]
        for rKmer, positions in d.items():
            score = ScoreKmers(qKmer, rKmer, matchScore, mismatchPen)
            if score >= threshHSSP:
                HSSPs[rKmer] = (positions, i)
    return HSSPs

# Input: a DNA string ref and an int k
# Output: a dictionary mapping each kmer in ref to a list of start positions 
def Indexation(ref, k):
    d = {}
    for i in range(len(ref)+1-k):
        kmer = ref[i:i+k]
        if kmer not in d:
            d[kmer] = [i]
        else:
            d[kmer].append(i)
    return d
"""

# Input: two DNA strings ref and query, an int k, scoring parameters and and HSSP threshold
# Output: a 2d list of seeds represented as [start coord in query, start coord in ref]
def BestSeeds(ref, query, k, matchScore, mismatchPen, matrix, threshHSSP):
    d = EncodedIndexation(ref, k)
    seeds = []
    allPossibleKmers = GenerateAllKmers(k)
    for i in range(len(query)+1-k):
        qKmer = query[i:i+k]

        # generate a list kmers, including qKmer, with an ungapped alignment at or above the HSSP threshold
        potentialHSSPs = []
        for kmer in allPossibleKmers:
            score = ScoreKmers(qKmer, kmer, matchScore, mismatchPen)
            if score >= threshHSSP:
                potentialHSSPs.append((kmer, score))

        # if a potential HSSP matches a reference kmer in d, save the qKmer start position (i) and the rKmer start position(s) (d[HSSP])
        for kmer, score in potentialHSSPs:
             check = KmerNumericalEncoding(kmer)
             if check in d:
                for pos in d[check]:
                    match = [i, pos, score]
                    seeds.append(match)          

        """
        # for each potential HSSP, search the bst for matches
        # for each match, save the qKmer start position (i) and the rKmer start position(s) (given in bst)
        for check in potentialHSSPs:
            # will return a tuple if a match is found
            match = SearchRefKmerBST(check, bst)
            if match != None:
                # might be multiple match positions in ref
                matchPositions = match[1]
                for pos in matchPositions:
                    coords = [i, pos]
                    seeds.append(coords)
        """

    return seeds

# Input: a reference sequence string and an int k
# Output: a dictionary mapping numerically-encoded kmers in ref to a list of their start positions in ref
def EncodedIndexation(ref, k):
    d = {}
    for i in range(len(ref)+1-k):
        kmer = ref[i:i+k]
        code = KmerNumericalEncoding(kmer)
        if code is None:
            continue
        if code not in d:
            d[code] = [i]
        else:
            d[code].append(i)
    return d

# Input: a kmer string
# Output: a numerical representation of the kmer
def KmerNumericalEncoding(kmer):
    encodingDict = {"A":0, "C":1, "G":2, "T":3}
    k = len(kmer)
    res = 0
    for pos in range(k):
        base = kmer[pos]
        # skip if not in encodingDict (like for N, other ambigious nucleotides)
        if base not in encodingDict:
            return None
        val = encodingDict[base]
        res += val * (4**(k-pos-1))
    return res

# Input: an int k
# Output: a list of all possible DNA strings of length k
def GenerateAllKmers(k):
    bases = ["A", "C", "T", "G"]
    return [''.join(p) for p in product(bases, repeat=k)]
    
# Input: two kmers of the same length and two scoring parameters
# Output: the score of the ungapped alignment between the kmers
def ScoreKmers(qKmer, rKmer, matchScore, mismatchPen):
    """
    Redudant with current scoreKmers function in extendSeeds.py, waiting for updated version that uses nucleotide BLOSUM equivalent
    """
    score = 0
    for i in range(len(qKmer)):
        if qKmer[i] == rKmer[i]:
            score += matchScore
        else:
            score -= mismatchPen
    return score 


"""
# Input: a dictionary mapping numerically-encoded kmers in ref to a list of their start positions in ref
# Output: a bst where nodes are key, value pairs in the dictionary
def ConstructRefKmerBST(d):
    tree = BST()
    for kmer in sorted(d):
        startPositions = d[kmer]
        if tree.root == None:
            tree.root = RefSeqNode(kmer)
            tree.root.pos = startPositions
        node = RefSeqNode
    pass
"""
