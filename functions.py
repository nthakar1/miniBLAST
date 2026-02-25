from bestSeeds import BestSeeds

# Parent function that returns a maximally scoring local alignment between 2 strings. CURRENTLY INCOMPLETE
# Input: (probably missing things currently)
# Output:
def miniBLAST(ref, query, k, matchScore, mismatchPen, threshHSSP):
    seeds = BestSeeds(ref, query, k, matchScore, mismatchPen, threshHSSP)
    # localAlignments = ExtendFromSeeds(ref, query, seeds, ...)
    pass