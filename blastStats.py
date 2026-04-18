import math
from datatypes import BlastConfig, BLASTN_PARAMS, BLASTP_PARAMS

def calculate_bit_score(S, params: BlastConfig) -> float:
    # Normalized score S' = (lambda * S - ln K) / ln 2 
    bit_score = (params.lambda_val * S - math.log(params.K)) / math.log(2)
    return bit_score

def calculate_e_value(S, m, n, params: BlastConfig) -> float:
    # the bit score is the normalized s' score
    # E = K * m * n * e^(-lambda * S)
    S_prime=calculate_bit_score(S,params)
    e_val = m * n * (2**-S_prime)
    
    return e_val

def compute_s1_threshold(query, db_records, params: BlastConfig, e_threshold=0.001):
    """
    Compute an ungapped-extension threshold s1 from the BLAST
    extreme-value approximation:

        E = K * m * n * exp(-lambda * S)

    Solving for S gives:

        s1 = ln(K * m * n / E_target) / lambda

    Parameters
    ----------
    query : str
        Query sequence.
    db_records : list
        List of SeqRecord objects in the database.
    params : BlastConfig
        Must contain:
            - params.K
            - params.lambda_val
    e_threshold : float, optional
        Target expected number of chance ungapped hits allowed through.
        Lower values make s1 stricter. Default = 0.001.

    Returns
    -------
    float
        Raw score threshold for deciding whether an ungapped hit
        should proceed to gapped extension.
    """
    if not query:
        raise ValueError("query must be a non-empty sequence")

    if not db_records:
        raise ValueError("db_records must be a non-empty list")

    if e_threshold <= 0:
        raise ValueError("e_threshold must be > 0")

    if params.K <= 0 or params.lambda_val <= 0:
        raise ValueError("params.K and params.lambda_val must both be > 0")

    m = len(query)
    n = sum(len(str(record.seq)) for record in db_records)

    if n <= 0:
        raise ValueError("database total length must be > 0")

    s1 = math.log((params.K * m * n) / e_threshold) / params.lambda_val
    return s1



