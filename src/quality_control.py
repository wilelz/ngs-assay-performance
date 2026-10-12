import numpy as np

def adaptive_trim(phred_scores, error_budget, min_length):

    error_probabilities = 10 ** (-phred_scores / 10)

    cumulative_errors = np.cumsum(error_probabilities)

    valid_positions = np.where(cumulative_errors <= error_budget)[0]

    if len(valid_positions) == 0:
        cutoff = 0
    else:
        cutoff = valid_positions[-1] + 1
        
    passed = cutoff >= min_length

    return int(cutoff), bool(passed)
