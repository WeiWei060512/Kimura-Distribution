# Kimura-Distribution

All mtDNA mutations initially affect a proportion of the many mtDNA molecules within a cell (heteroplasmy). Heteroplasmic pathogenic mtDNA mutations characteristically display a threshold effect, where single cells must harbour a high proportion of mutated molecules before they affect oxidative phosphorylation in the production of ATP. Individuals inheriting a high proportion of mtDNA heteroplasmy were more likely to have a severe disease than individuals with low levels of mtDNA heteroplasmy. 

Heteroplasmic mtDNA variants segregate rapidly during maternal transmission due a germ line genetic bottleneck, but the underlying mechanisms are largely unknown. Although random genetic drift is recognized as an important mechanism, selection mechanisms are thought to play a role as well.

To understand whether germline selection (besides random genetic drift) play a role during the transmission of heteroplasmic pathogenic mitochondrial DNA (mtDNA) mutations in humans, we tested a range of heteroplasmy levels in 96 quiescent oocytes from primordial follicles and 426 growing oocytes (from primary, secondary and antral follicles) isolated from 18 mice against the Kimura probability distribution [1]. 

Kimura probability distribution was firstly developed by a Japanese biologist Motoo Kimura for gene frequency probabilities in diploid populations. The theory has been used in mitochondrial heteroplasmy by the researchers.

This python module was inspired by Wonnapinij et al. [2]. It can be used to test the distribution of a data set of mtDNA heteroplasmy values for a group of individuals arising from a common founder against Kimura probability distribution.

The code includes the following steps:
Input: a text file contains a list of heteroplasmy values (between 0 to 1)

Step 1: Calculate mean, variance and b values (combined t (generations) and Neff (the effective population size)) from input heteroplasmy values (observed data).

Step 2: Fit the observed data into a Kimura probability distribution.

Step 3: Compare the distributions of observed data and derived Kimura probability distribution, and calculate D value using Kolmogorov-Smirnov (KS) test.

Step 4: Generate the random datasets using Monte-Carlo simulation.

Step 5: Calculate p value - the fraction of simulated data sets whose maximum deviation from the theoretical probability distribution was larger than the maximum deviation of the experimental data set.

[1] Selfish replication of mitochondrial DNA during oocyte maturation propagates mutation transmission (manuscript is under review).
[2] The distribution of mitochondrial DNA heteroplasmy due to random genetic drift. Wonnapinij P, Chinnery PF, Samuels DC. Am J Hum Genet. 2008 Nov;83(5):582-93.
