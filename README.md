# KAEstimate
Simulation-based parameter estimation for the Karlin-Altschul statistic (BLAST statistic).

## Description
The Karlin-Altschul statistic describes the probability to observe a local sequence alignment with a given score by chance in a database search. For this, it takes into account query and database length, letter frequencies and the scoring function. The distribution has two parameters that have to be computed based on the scoring function (mostly referred to as scoring matrix). For ungapped alignments, efficient numerical methods are described by Karlin and Altschul for computing said parameters. However, they cannot be directly applied to gapped alignments. This script computes the Karlin Altschul distribution parameters by simulating random local alignments.

The Karlin-Altschul formula gives the E-Value for an observed score S as: `E(S)=Kappa*m*n*exp((-lambda)*S)`. Where S corresponds to the score, m to the database sequence length, n to the query sequence length. The additional parameters lambda and Kappa can be computed using this script.

## Usage
Please see `example.py` for usage.

## Requirements
* Python 2.7+
* pairwise2 alignment package from BioPython
