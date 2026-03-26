from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds
from TwoHit import TwoHitSeeds

def miniBLAST(ref, query, k, matchReward, mismatchPen, gapOpenPen, gapExtendPen,
              threshHSSP, s1, Xdrop_ungap, Xdrop_gap, Xdrop_final, A, matrix=None):
    """
    Full miniBLAST pipeline:
      1. Seed generation  (BestSeeds)
      2. 2-Hit filter     (TwoHitSeeds)
      3. Extension        (ungapped → gapped with affine gaps)

    Parameters:
    ref, query     : DNA or protein strings
    k              : word / kmer length for seeding
    matchReward    : per-base match bonus (DNA; ignored when matrix is provided)
    mismatchPen    : per-base mismatch penalty (DNA; ignored when matrix is provided)
    gapOpenPen     : gap-open penalty
    gapExtendPen   : gap-extension penalty
    threshHSSP     : minimum kmer score to register as a seed
    s1             : minimum ungapped-extension score to pass into gapped phase
    Xdrop_ungap    : X-drop threshold during ungapped extension
    Xdrop_gap      : X-drop threshold during gapped extension
    Xdrop_final    : X-drop threshold applied during backtracking
    A              : maximum diagonal window for the 2-Hit filter
    matrix         : optional substitution matrix (e.g. BLOSUM62 for proteins)
    """
    singleSeeds = BestSeeds(ref, query, k, matchReward, mismatchPen, threshHSSP)
    print("Seeds1:", singleSeeds)

    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    print("Seeds2:", seedsToExtend)

    # Use 2-hit filtered seeds; fall back to all seeds if none survive the filter
    seeds_for_extension = seedsToExtend if seedsToExtend else singleSeeds

    bestLocalAlignment = extendFromSeeds(
        ref, query, seeds_for_extension, k,
        s1, matrix, matchReward, mismatchPen,
        gapOpenPen, gapExtendPen,
        Xdrop_ungap, Xdrop_gap, Xdrop_final
    )

    return bestLocalAlignment