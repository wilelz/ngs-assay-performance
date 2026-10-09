
from scipy.stats import binom


def calculate_sensitivity(depth, vaf, threshold):
    """Probability of detecting a true variant."""
    return binom.sf(threshold - 1, depth, vaf)


def calculate_specificity(depth, error_rate, threshold):
    """Probability of correctly rejecting a mutation-free sample."""
    return binom.cdf(threshold - 1, depth, error_rate)


def minimum_depth(
    vaf,
    error_rate,
    threshold,
    target_sensitivity=0.95,
    target_specificity=0.99,
    max_depth=10000
):       
    """Find the minimum depth meeting both performance targets."""  
    for depth in range(1, max_depth + 1):
        sensitivity = calculate_sensitivity(depth, vaf, threshold)
        specificity = calculate_specificity(depth, error_rate, threshold)

        if (sensitivity >= target_sensitivity
                and specificity >= target_specificity):
            return depth

    return None
 
def false_positive_probability(n_molecules, error_rate, threshold):
    """Find the probability of returning a false positive result."""
    probability = binom.sf(threshold - 1, n_molecules, error_rate)
    return probability
  

def molecular_detection_probability(n_molecules, vaf, threshold):
    """Find the probability of detecting a genuine mutant molecules."""
    probability = binom.sf(threshold-1, n_molecules, vaf)
    return probability
