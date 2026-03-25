from Bio.Align import substitution_matrices

# BLOSUM62 = industry standard protein scoring matrix, structured like a dictionary
BLOSUM62 = substitution_matrices.load("BLOSUM62")

# DNA-DNA alignment parameters (BLASTN)
BLASTN_PARAMS = {
    "name": "blastn",
    "matrix": None,  
    "match_reward": 2,
    "mismatch_penalty": 3,
    "gap_opening": 5,
    "gap_extension": 2,
    "k_mer_size": 7,  
    "xdrop_ungap": 20,
    "xdrop_gap": 30,
    "xdrop_gap_final": 100

}

# Protein-Protein alignment parameters (BLASTP)
BLASTP_PARAMS = {
    "name": "blastp",
    "matrix": BLOSUM62, 
    "gap_opening": 11, # protein gap penalties are higher
    "gap_extension": 1,
    "k_mer_size": 3,  # ~ codon length= 3
    "x_drop_ungapped": 7
}