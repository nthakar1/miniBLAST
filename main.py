from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds
from TwoHit import TwoHitSeeds

# Input:
# Output:
def miniBLAST(ref, query, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, threshHSSP, Xdrop, A):
    singleSeeds = BestSeeds(ref, query, k, matchReward, mismatchPen, threshHSSP)
    
    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    
    bestLocalAlignment = extendFromSeeds(ref, query, seedsToExtend, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, Xdrop)
    
    return bestLocalAlignment