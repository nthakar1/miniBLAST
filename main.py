from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds
from TwoHit import TwoHitSeeds
from datatypes import BLASTN_PARAMS, BLASTP_PARAMS, BlastConfig

# Input:
# Output:
def miniBLASTn(ref, query, threshHSSP, A):
    k=BLASTN_PARAMS.k_mer_size
    matchReward=BLASTN_PARAMS.match_reward
    mismatchPen=BLASTN_PARAMS.mismatch_penalty
    gapOpenPen=BLASTN_PARAMS.gap_opening
    gapExtendPen=BLASTN_PARAMS.gap_extension
    Xdrop=BLASTN_PARAMS.xdrop_gap_final

    singleSeeds = BestSeeds(ref, query, k, matchReward, mismatchPen, threshHSSP)
    print("Seeds1:", singleSeeds)
    
    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    print("Seeds2:", seedsToExtend)

    bestLocalAlignment = extendFromSeeds(ref, query, singleSeeds, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen, Xdrop)
    
    return bestLocalAlignment