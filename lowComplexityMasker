"""lowComplexityMasker used in BLAST+ to mask repetitive regions. 
Using this in miniBLAST because viral genomes tend to have repetitive regions, 
making alignment algorithm much slower."""


try:
    import pydustmasker
    PYDUSTMASKER_AVAILABLE = True
except ImportError:
    PYDUSTMASKER_AVAILABLE = False
    print("Warning: pydustmasker not available. Install with: pip install pydustmasker")

def mask_low_complexity_regions(sequence, score_threshold=20, mask_char='lowercase', window_size=64):
    """ Uses NCBI standard symmetrical dust (SDUST) to mask low complexity regions with repeats. 
    Repetitive regions are masked by converting the uppercase nucleotide sequence to a lowercase sequence.
    -----
    Input: sequence (str)-- DNA sequence that may be masked
    -----
    Parameters: score_threshold (int)-- default= 20, any regions scoring above 20 are classified as a low complexity region and get masked
    window length (int)-- default= 64, window over which we should assess frequency of repeats
    mask_char : str, default 'lowercase'
    -----
    Output: sequence with masked regions (nucleotides are lowercase) and their intervals
    """
    if not PYDUSTMASKER_AVAILABLE:
        return sequence, []
    
    if not sequence or not isinstance(sequence, str):
        return sequence, []
        
    try:
        # Use SDUST algorithm (same as NCBI BLAST)
        masker = pydustmasker.DustMasker(
            sequence, 
            score_threshold=score_threshold,
            window_size=window_size
        )
        
        intervals = list(masker.intervals)
        
        # if the value for mask_char input was not in lowercase already, the lower method makes it lowercase
        if mask_char.lower() == 'lowercase':
            # soft mask method will convert everything to lowercase (soft mask by default when mask(hard: bool = True) is not set)
            masked_sequence = masker.mask()
        else:
            # Alternative: replace with N's
            masked_sequence = masker.mask(hard=True)
        return masked_sequence, intervals
        
    except Exception as e:
        print(f"Error in low complexity masking: {e}")
        return sequence, []

def masking_stats(sequence, intervals):
    """
    Get statistics about masking results.
    
    Returns:
    --------
    dict with keys: total_length, masked_bases, pct_masked, n_intervals
    """
    if not intervals:
        return {
            'total_length': len(sequence),
            'masked_bases': 0,
            'pct_masked': 0.0,
            'n_intervals': 0
        }
    
    # sum the lengths of every interval to determine the number of masked bases
    masked_bases = sum(end - start for start, end in intervals) 
    
    return {
        'total_length': len(sequence),
        'masked_bases': masked_bases, # sum of masked bases
        'pct_masked': round(100.0 * masked_bases / len(sequence), 2), # how many bases were masked in the entire sequence
        'n_intervals': len(intervals) # how many intervals are there
    }

def filter_seeds_by_masking(seeds, k, query_masked_intervals, ref_masked_intervals, overlap_threshold=0.5):
    """
    Filter seeds that significantly overlap with masked regions.
    -----------
    Inputs:
    seeds (list): List of [query_start, ref_start, score] seeds
    k (int): Kmer length  
    query_masked_intervals, ref_masked_intervals (list of tuples): Masked intervals
    overlap_threshold (float): default=0.5-- Fraction of kmer that can overlap with masked regions

    -------- 
    Output:
    list of filtered_seeds with overlap length<= k*overlap_threshold
    """

    # if there aren't low complexity regions in the query nor the ref, return the original seeds
    if not query_masked_intervals and not ref_masked_intervals:
        return seeds
    
    # number of masked bases per window = k*overlap_threshold. Example: k=3, overlap_threshold=0.5, max number of masked bases per window = 1.5
    max_masked_bases = int(k * overlap_threshold)
    filtered = []
    
    # for every seed in the list
    for query_start, ref_start, score in seeds:
        query_end = query_start + k
        ref_end = ref_start + k
        
        masked_bases = 0
        
        # determine indices for where the overlap occurs between the query seeds and the query masked intervals
        for start, end in query_masked_intervals:
            overlap_start = max(query_start, start)
            overlap_end = min(query_end, end)
            # determine whether the start of the overlap region does in fact occur at an index before the end of the overlap region (otherwise it isn't an overlap region)
            if overlap_start < overlap_end:
                # calculate number of nucleotides in overlap region
                masked_bases += overlap_end - overlap_start
                
       # determine indices for where the overlap occurs between the ref seeds and the ref masked intervals
        for start, end in ref_masked_intervals:
            overlap_start = max(ref_start, start)
            overlap_end = min(ref_end, end)
            if overlap_start < overlap_end:
                masked_bases += overlap_end - overlap_start

        # add a seed into our list of seeds only if it is less than the number of allowed overlapping bases between the intervals and the seeds
        if masked_bases <= max_masked_bases:
            filtered.append([query_start, ref_start, score])
            
    return filtered


def demo_masking():
    """Demonstrate DUST masking on test sequences."""
    
    if not PYDUSTMASKER_AVAILABLE:
        print("pydustmasker not available for demo")
        return
        
    # Test sequences with different complexity levels
    test_sequences = [
        ("High complexity", "ATGCAGTTCGAACCTGATCCGAATGCTA"),
        ("Poly-A (very low complexity)", "AAAAAAAAAAAAAAAAAAAA"),
        ("Dinucleotide repeat", "ATCGATATATATATATATATCGATCG"),
        ("Mixed complexity", "ATCGATCGAAAAAAAACGTACGTACGTTTTTTTTATCG")
    ]
    
    print("DUST Masking Demo")
    print("=" * 50)
    
    for name, seq in test_sequences:
        print(f"\n{name}:")
        print(f"Original:  {seq}")
        
        # Test with default threshold (20)
        masked, intervals = mask_low_complexity_regions(seq, score_threshold=20)
        stats = masking_stats(seq, intervals)
        print(f"Masked:    {masked}")
        print(f"Intervals: {intervals}")
        print(f"Masked:    {stats['pct_masked']:.1f}% ({stats['masked_bases']}/{stats['total_length']} bp)")

if __name__ == "__main__":
    demo_masking()



