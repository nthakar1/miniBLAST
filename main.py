from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds
from TwoHit import TwoHitSeeds
from datatypes import BLASTN_PARAMS, BLASTP_PARAMS, BlastConfig
from database import fetch_mixed_database

from Bio import SeqIO
import pandas as pd
import time

def miniBLASTn(ref, query, s1, A):
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
    k=BLASTN_PARAMS.k_mer_size
    matchReward=BLASTN_PARAMS.match_reward
    mismatchPen=BLASTN_PARAMS.mismatch_penalty
    matrix=BLASTN_PARAMS.matrix
    gapOpenPen=BLASTN_PARAMS.gap_opening
    gapExtendPen=BLASTN_PARAMS.gap_extension
    Xdrop_ungap=BLASTN_PARAMS.xdrop_ungap
    Xdrop_gap=BLASTN_PARAMS.xdrop_gap
    Xdrop_final=BLASTN_PARAMS.xdrop_gap_final

    seedMismatchAllowed = int(k/3) # should be in the range of 1-3 mismatches for the expected inputs of k
    threshHSSP = matchReward*(k-seedMismatchAllowed) + mismatchPen*(seedMismatchAllowed)

    startTime = time.perf_counter()
    singleSeeds = BestSeeds(ref, query, k, matchReward, mismatchPen, threshHSSP)
    #print("Seeds1:", singleSeeds)
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print("Seed identification complete. Time elapsed:", (endTime-startTime)/60, "seconds")

    startTime = time.perf_counter()
    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    #print("Seeds2:", seedsToExtend)
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print("Two Hit filtration complete. Time elapsed:", (endTime-startTime)/60, "seconds")

    # Use 2-hit filtered seeds; fall back to all seeds if none survive the filter
    seeds_for_extension = seedsToExtend if seedsToExtend else singleSeeds
    
    startTime = time.perf_counter()
    bestLocalAlignment = extendFromSeeds(
        ref, query, seeds_for_extension, k,
        s1, matrix, matchReward, mismatchPen,
        gapOpenPen, gapExtendPen,
        Xdrop_ungap, Xdrop_gap, Xdrop_final
    )
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print("Local alighment complete. Time elapsed:", (endTime-startTime)/60, "seconds")

    return bestLocalAlignment

def main():
    
    queriesViral = [
        # H1N1 sequences (your target of interest)
        {
            "query": "H1N1[All Fields] AND influenza A virus[Organism]",
            "max_results": 200,
            "label": "H1N1_HA"
        },
        {
            "query": "H1N1[All Fields] AND influenza A virus[Organism]", 
            "max_results": 200,
            "label": "H1N1_NA"
        },
        
        # Other influenza subtypes (closely related background)
        {
            "query": "H3N2[All Fields] AND influenza A virus[Organism]",
            "max_results": 150,
            "label": "H3N2_HA"
        },
        {
            "query": "H5N1[All Fields] AND influenza A virus[Organism]",
            "max_results": 100,
            "label": "H5N1"
        },
        
        # Other respiratory viruses (distant background for E-value diversity)
        {
            "query": "SARS-CoV-2[Organism] AND complete genome[Title]",
            "max_results": 100,
            "label": "COVID"
        },
        {
            "query": "rhinovirus[Organism] AND complete genome[Title]",
            "max_results": 100,
            "label": "rhinovirus"
        },
        
        # Some bacterial background (adds statistical noise)
        {
            "query": "16S ribosomal RNA[Gene] AND bacteria[Organism]",
            "max_results": 150,
            "label": "bacteria_16S"
        }
    ]
    
    segments = []
    query_file = "human_h1n1_1934.fna"
    for record in SeqIO.parse(query_file, "fasta"):
        segments.append(str(record.seq))
    print(f"Found {len(segments)} sequence segments in query file")
    
    # could consider computing once and storing--time vs memory
    db_h1n1 = fetch_mixed_database(queriesViral, "h1n1_mixed_database.fasta")

    query = segments[0]
    alignments = []

    # for test: only aligning to 10 of the roughly 1000 sequences in the database
    for i in range(0, len(db_h1n1), 100):
        print(f"Aligning to reference {i}")

        record = db_h1n1[i]
        ref = str(record.seq)

        # each alignment is the best alignment between the query and that ref, as a dictionary containing score, alignment, position, and coverage
        ### later: is position position in reference? where are % identity and stat values being output
        alignment = miniBLASTn(query, ref, 0, float("inf"))
        alignments.append(alignment)

    df = pd.DataFrame(alignments)
    df.to_csv("blast_results.csv")

if __name__ == "__main__":
    main()