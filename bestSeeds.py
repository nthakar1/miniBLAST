from itertools import product
from extendSeeds import scorePair
import time

##### NOTES FOR LATER: #####
"""
- Will want to account for IUPAC nucleotide ambiguity codes during indexation and scoring (reference or query seqs with degenerated bases)
"""

# Input: two DNA strings ref and query, an int k, scoring parameters and and HSSP threshold
# Output: a 2d list of seeds represented as [start coord in query, start coord in ref]
def BestSeeds(ref, query, k, matchScore, mismatchPen, matrix, threshHSSP):
    d = EncodedIndexation(ref, k)
    seeds = []
    allPossibleKmers = GenerateAllKmers(k) # move outside this function?
    for i in range(len(query)+1-k):
        qKmer = query[i:i+k]

        # generate a list kmers, including qKmer, with an ungapped alignment at or above the HSSP threshold
        potentialHSSPs = [] # FIND A DIFFERENT WAY TO GENERATE THIS
        for kmer in allPossibleKmers:
            score = ScoreKmers(qKmer, kmer, matrix, matchScore, mismatchPen)
            if score >= threshHSSP:
                potentialHSSPs.append((kmer, score))

        # if a potential HSSP matches a reference kmer in d, save the qKmer start position (i) and the rKmer start position(s) (d[HSSP])
        for kmer, score in potentialHSSPs:
             check = KmerNumericalEncoding(kmer)
             if check in d:
                for pos in d[check]:
                    match = [i, pos, score]
                    seeds.append(match)          

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
def ScoreKmers(qKmer, rKmer, matrix, matchScore, mismatchPen):
    score = 0
    for i in range(len(qKmer)):
        score += scorePair(qKmer[i], rKmer[i], matrix, matchScore, mismatchPen)
        if qKmer[i] == rKmer[i]:
            score += matchScore
        else:
            score -= mismatchPen

    return score 

