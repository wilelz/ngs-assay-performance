# NGS Assay Performance Analysis

## Project Overview

This project investigates the statistical limits of detecting low-frequency genetic variants using next-generation sequencing (NGS).

Using Python, binomial probability models, and computational simulations, the project examines how sequencing depth, variant allele frequency (VAF), molecular sampling, background error rates, and detection thresholds affect theoretical assay sensitivity and specificity.

The goal is to develop reproducible scientific computing workflows for evaluating rare-variant detection strategies.

**Project status:** Exploratory computational modeling with reusable Python functions, automated tests, and GitHub Actions continuous integration.

## Research Questions

1. How do sequencing depth and VAF influence the probability of detecting a mutation?
2. How do sequencing errors affect sensitivity and false-positive rates?
3. Why does unique molecular depth differ from raw sequencing depth?
4. How do detection thresholds change the trade-off between sensitivity and specificity?
5. What minimum number of independent DNA molecules is required to satisfy specified performance targets?

## Methods

### 1. Read-Level Variant Detection

Variant-supporting read counts are modeled using a binomial distribution.

For sequencing depth \(n\), variant allele frequency \(p\), and detection threshold \(k\):

\[
P(X \geq k) = P(\operatorname{Binomial}(n,p) \geq k)
\]

The model evaluates detection probability and specificity across different depths, allele frequencies, and thresholds.

### 2. Sequencing Error Modeling

The analysis investigates background sequencing errors and their effects on apparent variant-supporting reads.

An exploratory error-aware model accounts for both true variant observations and erroneous reference-to-alternate observations under simplified symmetric error assumptions.

### 3. Molecular Sampling

The project distinguishes between total sequencing reads and independent original DNA molecules.

The probability of sampling at least one mutant molecule is:

\[
P(X \geq 1)=1-(1-\mathrm{VAF})^N
\]

where \(N\) is the number of independent DNA molecules.

The model also calculates the probability of sampling at least a specified number of mutant molecules.

### 4. Molecular Detection Thresholds

The molecular model calculates:

- Sensitivity: probability that a mutation-positive sample meets the molecular detection threshold.
- False-positive probability: probability that background molecular errors meet the threshold in a mutation-free sample.
- Specificity: one minus the false-positive probability.

These calculations use independent binomial sampling assumptions.

### 5. Minimum Molecular Sample Size

A parameterized Python function searches across candidate molecular sample sizes and detection thresholds to identify the first combination satisfying specified sensitivity and specificity targets.

The current search evaluates molecule counts from 3,000 to 20,000 and thresholds from 1 to 10.

## Key Results

For a hypothetical assay with:

| Parameter | Value |
|---|---|
| True VAF | 0.1% |
| Molecular error probability | 0.01% |
| Minimum sensitivity | 90% |
| Minimum specificity | 99% |
| Minimum qualifying sample size | **6,679 molecules** |
| Qualifying detection threshold | **4 independent mutant molecules** |

The search identified 6,679 as the first qualifying
