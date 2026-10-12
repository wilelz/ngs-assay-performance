import numpy as np
from src.quality_control import adaptive_trim

def test_adaptive_trim_high_quality():
    phred_scores = np.full(200, 40)

    cutoff, passed = adaptive_trim(phred_scores, 1.0, 150)
    assert cutoff == 200
    assert passed 

def test_adaptive_trim_low_quality():
    phred_scores = np.full(200, 10)

    # Call adaptive_trim()
    cutoff, passed = adaptive_trim(phred_scores, 1.0, 150)
    assert cutoff == 10
    assert not passed 

def test_adaptive_trim_150_quality():
    phred_scores = np.full(150, 30)

    # Call adaptive_trim()
    cutoff, passed = adaptive_trim(phred_scores, 1.0, 150)
    assert cutoff == 150
    assert passed
