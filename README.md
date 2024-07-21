# Double Gene Citations as Randomizing Factor When Constructing Genetic Bayesian Networks

This is a script I wrote during my internship in early 2021

The Goal: The number of times a gene pair is found to be associated together in the same research paper's abstract
is used as a factor when randomally picking gene pairs as priors to train a genetic Bayesian Network.  

1. genes expressing in a particula organ are taken from the Human Protein Atlas DB.
2. Each gene's single citation number, the number of times the gene was found in research paper's abstracts, is taken from Pubmed.
3. Each gene pair's double citation number - the number of times the two genes were found in research paper's abstracts together, is taken from Pubmed.
4. Relative associatin factor is calculated.
5. Results are binned into ten bins which control the chance a gene pair wil be randomally picked to construct the naschent bayesian network during its training.
