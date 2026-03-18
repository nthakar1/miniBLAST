from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds
from TwoHit import TwoHitSeeds

# Input:
# Output:
def miniBLAST(ref, query, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, threshHSSP, Xdrop, A):
    singleSeeds = BestSeeds(ref, query, k, matchReward, mismatchPen, threshHSSP)
    print("Seeds1:", singleSeeds)
    
    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    print("Seeds2:", seedsToExtend)

    bestLocalAlignment = extendFromSeeds(ref, query, singleSeeds, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, Xdrop)
    
    return bestLocalAlignment