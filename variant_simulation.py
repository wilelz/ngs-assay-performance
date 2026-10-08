import numpy as np

# Parameters
sequencing_depth = 500
variant_frequency = 0.01

# Simulate one sequencing experiment
variant_reads = np.random.binomial(
    n=sequencing_depth,
    p=variant_frequency
)

print(f"Sequencing depth: {sequencing_depth}")
print(f"Variant-supporting reads: {variant_reads}")
