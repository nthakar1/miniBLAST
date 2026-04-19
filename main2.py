from bestSeeds import BestSeeds
from extendSeeds import extendFromSeeds
from TwoHit import TwoHitSeeds
from datatypes import BLASTN_PARAMS, BLASTP_PARAMS, BlastConfig, BLOSUM
from database import fetch_mixed_database
from blastStats import calculate_bit_score, calculate_e_value, compute_s1_threshold

from Bio import SeqIO
import csv
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

    seedMismatchAllowed = 1
    threshHSSP = matchReward*(k-seedMismatchAllowed) - mismatchPen*(seedMismatchAllowed)

    startTime = time.perf_counter()
    singleSeeds, masking_info = BestSeedsWithMasking(
        ref=ref, query=query, k=k, matchScore=matchReward, mismatchPen=mismatchPen, matrix=None, threshHSSP=threshHSSP,
        apply_masking=True, dust_complexity_threshold=20
    )
    
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print(f"{len(singleSeeds)} seeds identified")
    if len(singleSeeds) == 0:
        print("No seeds identified in intial search. Haulting run.")
        return
    print("Seed identification complete. Time elapsed:", round((endTime-startTime), 5), "seconds")

    startTime = time.perf_counter()
    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    #print("Seeds2:", seedsToExtend)
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print(f"{len(seedsToExtend)} seeds remaining after filtration")
    if len(seedsToExtend) == 0:
        print("No seeds left to extend after two-hit filtration. Haulting run.")
        return
    print("Two Hit filtration complete. Time elapsed:", round((endTime-startTime), 5), "seconds")

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
    print("Local alighment complete. Time elapsed:", round((endTime-startTime), 5), "seconds")

    return bestLocalAlignment

def miniBLASTp(ref, query, s1, A):
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
    k=BLASTP_PARAMS.k_mer_size
    matchReward=BLASTP_PARAMS.match_reward
    mismatchPen=BLASTP_PARAMS.mismatch_penalty
    matrix=BLASTP_PARAMS.matrix
    gapOpenPen=BLASTP_PARAMS.gap_opening
    gapExtendPen=BLASTP_PARAMS.gap_extension
    Xdrop_ungap=BLASTP_PARAMS.xdrop_ungap
    Xdrop_gap=BLASTP_PARAMS.xdrop_gap
    Xdrop_final=BLASTP_PARAMS.xdrop_gap_final
    threshT=BLASTP_PARAMS.threshold_T

    startTime = time.perf_counter()
    singleSeeds = BestSeeds(ref, query, k, matchReward, mismatchPen, matrix, threshT)
    
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print(f"{len(singleSeeds)} seeds identified")
    if len(singleSeeds) == 0:
        print("No seeds identified in intial search. Haulting run.")
        return
    print("Seed identification complete. Time elapsed:", round((endTime-startTime), 5), "seconds")

    startTime = time.perf_counter()
    seedsToExtend = TwoHitSeeds(singleSeeds, k, A)
    #print("Seeds2:", seedsToExtend)
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print(f"{len(seedsToExtend)} seeds remaining after filtration")
    if len(seedsToExtend) == 0:
        print("No seeds left to extend after two-hit filtration. Haulting run.")
        return
    print("Two Hit filtration complete. Time elapsed:", round((endTime-startTime), 5), "seconds")

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
    print("Local alighment complete. Time elapsed:", round((endTime-startTime), 5), "seconds")

    return bestLocalAlignment

def main():
    
    startTime = time.perf_counter()
    
    # 1. CHANGED: Updated the queries to the specific HA database search
    queriesViral = [
        {
            "query": 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H1N1[All Fields]',
            "max_results": 25,
            "label": "H1N1_HA"
        },
        {
            "query": 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H3N2[All Fields]',
            "max_results": 25,
            "label": "H3N2_HA"
        },
        {
            "query": 'Influenza A virus[Organism] AND ("hemagglutinin"[Title] OR "segment 4"[Title]) AND H5N1[All Fields]',
            "max_results": 25,
            "label": "H5N1_HA"
        },
        {
            "query": 'Influenza B virus[Organism] AND hemagglutinin[Title]',
            "max_results": 15,
            "label": "InfluenzaB_HA"
        }
    ]

    segments = []
    # Set the input sequence to local single HA file
    query_file = "H1N1_HA_single.fasta" 
    
    for record in SeqIO.parse(query_file, "fasta"):
        segments.append(str(record.seq))
    print(f"Found {len(segments)} sequence segments in query file")
    
    # could consider computing once and storing--time vs memory
    db_h1n1 = fetch_mixed_database(queriesViral, "nucleotide", "viral_HA_mixed_database.fasta")

    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print()
    print("Imported database and query sequence. Time elapsed:", round((endTime-startTime), 5), "seconds")

    query = segments[0]
    s1 = compute_s1_threshold(query, db_h1n1, BLASTN_PARAMS)
    print(f"Computed s1 threshold: {round(s1, 2)}")
    alignments = []
    ref_descriptions = []

    # for test: only aligning to 10 of the roughly 1000 sequences in the database
    startTime = time.perf_counter()
    
    # NOTE: Since your new queries list only pulls ~90 total sequences, a step 
    # of 100 here means it will only run on the first sequence (i=0). 
    # If you want it to test more, change the 100 to a 10 or a 1!
    for i in range(0, len(db_h1n1), 10):
        print()
        print(f"Aligning to reference {i}")

        record = db_h1n1[i]
        ref = str(record.seq)

        source = record.description
        ref_descriptions.append(source)

        # each alignment is the best alignment between the query and that ref, as a dictionary containing score, alignment, position, and coverage
        alignment = miniBLASTn(ref, query, s1=s1, A=40)
        
        if alignment:
            # ---> ADDED: Store the reference description in the dictionary <---
            alignment["reference_sequence"] = source
            
            alignment["bit_score"] = round(calculate_bit_score(alignment["score"], BLASTN_PARAMS), 4)
            alignment["e_value"]   = round(calculate_e_value(alignment["score"], len(query), len(db_h1n1), BLASTN_PARAMS), 4)
            
            alignments.append(alignment)

    
    sum(i**2 for i in range(1000000))
    endTime = time.perf_counter()
    print()
    print("Finished aligning to database. Time elapsed:", round((endTime-startTime), 5)/60, "minutes")

    alignments.sort(key=lambda x: x["bit_score"], reverse=True)

    # ---> CHANGED: Added "reference_sequence" as the first fieldname <---
    with open("blast_HA.csv", "w", newline="") as f:
      w = csv.DictWriter(f, fieldnames=["reference_sequence", "score", "bit_score", "e_value", "alignment", "position", "query_coverage", "pct_identity"])
      w.writeheader() 
      w.writerows(a for a in alignments if a is not None)
    print("Results written to csv!")

if __name__ == "__main__":
    main()
