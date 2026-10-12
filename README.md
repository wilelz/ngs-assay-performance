# NGS Assay Performance Analysis

## Project Overview

This project investigates the statistical limits of detecting low-frequency genetic variants using next-generation sequencing (NGS).

Using Python, binomial probability models, and computational simulations, the project examines how sequencing depth, variant allele frequency (VAF), molecular sampling, background error rates, and detection thresholds affect theoretical assay sensitivity and specificity.

The goal is to develop reproducible scientific computing workflows for evaluating rare-variant detection strategies.

**Project status:** Exploratory theoretical modeling (Version 0.1) and real FASTQ quality-control analysis (Version 0.2), with reusable Python functions, automated tests, and GitHub Actions continuous integration.

## Version 0.1 — Theoretical Assay Performance Modeling

### Research Questions

1. How do sequencing depth and VAF influence the probability of detecting a mutation?
2. How do sequencing errors affect sensitivity and false-positive rates?
3. Why does unique molecular depth differ from raw sequencing depth?
4. How do detection thresholds change the trade-off between sensitivity and specificity?
5. What minimum number of independent DNA molecules is required to satisfy specified performance targets?

### Methods

#### 1. Read-Level Variant Detection

Variant-supporting read counts are modeled using a binomial distribution.

For sequencing depth \(n\), variant allele frequency \(p\), and detection threshold \(k\):

\[
P(X \geq k) = P(\operatorname{Binomial}(n,p) \geq k)
\]

The model evaluates detection probability and specificity across different depths, allele frequencies, and thresholds.

#### 2. Sequencing Error Modeling

The analysis investigates background sequencing errors and their effects on apparent variant-supporting reads.

An exploratory error-aware model accounts for both true variant observations and erroneous reference-to-alternate observations under simplified symmetric error assumptions.

#### 3. Molecular Sampling

The project distinguishes between total sequencing reads and independent original DNA molecules.

The probability of sampling at least one mutant molecule is:

\[
P(X \geq 1)=1-(1-\mathrm{VAF})^N
\]

where \(N\) is the number of independent DNA molecules.

The model also calculates the probability of sampling at least a specified number of mutant molecules.

#### 4. Molecular Detection Thresholds

The molecular model calculates:

- Sensitivity: probability that a mutation-positive sample meets the molecular detection threshold.
- False-positive probability: probability that background molecular errors meet the threshold in a mutation-free sample.
- Specificity: one minus the false-positive probability.

These calculations use independent binomial sampling assumptions.

#### 5. Minimum Molecular Sample Size

A parameterized Python function searches across candidate molecular sample sizes and detection thresholds to identify the first combination satisfying specified sensitivity and specificity targets.

The current search evaluates molecule counts from 3,000 to 20,000 and thresholds from 1 to 10.

### Key Results

For a hypothetical assay with:

| Parameter | Value |
|---|---|
| True VAF | 0.1% |
| Molecular error probability | 0.01% |
| Minimum sensitivity | 90% |
| Minimum specificity | 99% |
| Minimum qualifying sample size | **6,679 molecules** |
| Qualifying detection threshold | **4 independent mutant molecules** |

The search identified **6,679 independent molecules** as the first qualifying sample size within the evaluated search range, with a detection threshold of four mutant molecules. This is a hypothetical binomial-model result, not a validated clinical assay.

**Notebook:** [01_variant_detection.ipynb](01_variant_detection.ipynb)

## Version 0.2 — Real FASTQ Quality Analysis and Adaptive Trimming

### Objective and Dataset

We analyzed **1,000 real Illumina MiSeq reads** (301 bases each) from the first FASTQ mate of publicly available ENA run **ERR4387389**, using a partial download rather than the complete sequencing run. The objective was to characterize read-position quality, estimate base-call error probabilities from Phred scores, and compare fixed-length and adaptive trimming.

**Notebook:** [02_real_fastq_analysis.ipynb](02_real_fastq_analysis.ipynb)

### Methods

- Decoded FASTQ Phred+33 quality scores using `ord(character) - 33`.
- Converted each score Q to its estimated base-error probability using `10 ** (-Q / 10)`.
- Characterized position-dependent mean Phred quality and proportions of bases meeting Q20 and Q30.
- Compared fixed read-prefix trimming with per-read adaptive trimming. Adaptive trimming retains the longest prefix whose cumulative estimated base-error probability is within a user-specified budget, then rejects qualifying prefixes below a minimum acceptable length.
- Implemented the reusable `adaptive_trim()` function and automated tests for high-quality, low-quality, and minimum-length inputs.

The reported adaptive analysis uses an expected-error budget of **1.0 per read** and a minimum acceptable retained length of **150 bases**.

### Results

| Measure | Result |
|---|---:|
| Original reads / bases | 1,000 / 301,000 |
| Reads passing adaptive length and error constraints | 943 (94.3%) |
| Bases retained by adaptive trimming | 231,361 (76.9%) |
| Mean length of surviving reads | 245.3 bases |
| Original total Phred-estimated errors | 5,267.6 |
| Adaptive total Phred-estimated errors | 755.9 |
| Reduction in total estimated errors versus original | **85.7%** |

#### Fixed 210-base versus adaptive trimming

| Measure | Fixed 210 bases | Adaptive, minimum 150 bases |
|---|---:|---:|
| Reads retained | 1,000 | 943 |
| Bases retained | 210,000 | 231,361 |
| Total estimated base-call errors | 988.41 | 755.9 |
| Mean estimated errors per retained read | 0.988 | 0.802 |

In this sample, adaptive trimming retained **10.17% more bases** and yielded **23.52% fewer total Phred-estimated errors** compared with fixed trimming at 210 bases. This is not a fully matched comparison because the adaptive strategy discards 57 short reads.

### Limitations

Phred-derived probabilities are estimated base-call errors, not experimentally verified mismatches or false-positive variant calls. A budget of 1.0 constrains the **expected** errors, not the realized errors. The chosen budget and 150-base minimum are exploratory; no reference alignment, locus-specific coverage assessment, variant calling, or paired-end analysis has yet validated these trimming decisions. The sampled reads may not represent the entire sequencing run, and read-level errors are not interchangeable with independent molecular error rates used in Version 0.1.

## Repository Structure

| File or directory | Purpose |
|---|---|
| `01_variant_detection.ipynb` | Theoretical rare-variant detection analysis |
| `02_real_fastq_analysis.ipynb` | Real FASTQ quality and trimming analysis |
| `src/assay_performance.py` | Theoretical model functions |
| `src/quality_control.py` | Adaptive-trimming function |
| `tests/` | Automated Python tests |
| `requirements.txt` | Python dependencies |
| `.github/workflows/tests.yml` | Continuous-integration tests |

## Installation and Testing

From the repository root:

```bash
python -m pip install -r requirements.txt
python -m pytest
```

## Future Work

Evaluate reference alignment, paired-end read handling, usable variant-site coverage, and rare-variant calling performance under different trimming thresholds.
