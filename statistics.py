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




