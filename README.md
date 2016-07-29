# KAEstimate
Simulation-based parameter estimation for the Karlin-Altschul distribution (BLAST statistic).

##Description
The Karlin-Altschul distribution gives the probability to observe a local sequence alignment with a given score by chance in a database search. For this, it takes into account query and database length, letter frequencies and the scoring function. The distribution has two parameters that have to be computed based on the scoring function. This script accomplishes this by simulating random local alignments and successively fititng the Karlin-Altschul parameters so that the distribution fits the observed score frequencies.

The Karlin-Altschul formula gives the E-Value for an observed score S as: `p(S)=Kappa*m*n*exp((-lambda)*S)`. Where S corresponds to the score, m to the database sequence length, n to the query sequence length. The additional parameters lambda and Kappa can be computed using this script.

##Usage
Please see `example.py` for usage.

##Requirements
* Python 2.7+
* pairwise2 alignment package from BioPython