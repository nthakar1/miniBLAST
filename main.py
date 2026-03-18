from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds

# Input:
# Output:
def miniBLAST(ref, query, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, threshHSSP, Xdrop):
    seeds = BestSeeds(ref, query, k, matchReward, mismatchPen, threshHSSP)
    bestLocalAlignment = extendFromSeeds(ref, query, seeds, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, Xdrop)
    pass