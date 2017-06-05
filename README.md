# DistMCMC
This is my coursework of Molecular Evolution and Phylogenetics (by guest lecturer Dr. Xuhua Xia from University of Ottawa, Canda, in summer 2014), which received the highest score among a class of about twenty students.
Estimation of the number of nucleotide substitution between DNA sequences is a fundamental question in the field of molecular evolution. The JC69 and K80 are the simplest nucleotide substitution models. Maximum likelihood estimation is a method of estimating parameters under a certain model. Markov chain Monte Carlo can be used to get the probability distribution of parameters by sampling from the equilibrium distribution of a Markov chain that has experienced enough generations. To learn nucleotide substitution models, maximum likelihood estimation, and Markov chain Monte Carlo, I write a Perl program called DistMCMC, which perform maximum likelihood estimation of distance between two DNA sequences using Markov chain Monte Carlo. For learning purpose, only JC69, JC69+Γ, K80 and K80+Γ are implemented.

On May 6, 2017, I realized that I had misunderstood some basic concepts in the coursework. The perl script (Yang2006exercise54.pl) of the original exercise in Yang 2006 was added to this folder.
