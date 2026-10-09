
import sys
from pathlib import Path

# Import our analysis module from src/
sys.path.insert(
    0,
    str(Path(__file__).resolve().parents[1] / "src")
)

from assay_performance import (
    calculate_sensitivity,
    calculate_specificity,
    minimum_depth,
    false_positive_probability,
    molecular_detection_probability
)


def test_sensitivity_increases_with_depth():
    """More coverage should improve sensitivity."""
    low = calculate_sensitivity(100, 0.01, 3)
    high = calculate_sensitivity(500, 0.01, 3)

    assert high > low


def test_specificity_increases_with_threshold():
    """Requiring more supporting reads reduces false positives."""
    low = calculate_specificity(500, 0.001, 1)
    high = calculate_specificity(500, 0.001, 3)

    assert high > low


def test_minimum_depth():
    """Check our previously calculated result."""
    result = minimum_depth(
        vaf=0.01,
        error_rate=0.0005,
        threshold=3
    )

    assert result == 628


def test_molecular_detection_probability():
    result = molecular_detection_probability(3000, 0.001, 3)

    assert round(result, 4) == 0.5769

