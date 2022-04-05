# Coala model 

# usage coal_model(sample_size, loci_number = 0, loci_length = 1000, ploidy = 1)
#coal_model(76, loci_number = 1,ploidy = 2)


# A model with one population and 1 loci:
#model <- coal_model(sample_size = 76, loci_number = 1, ploidy = 2) +
 # feat_mutation(4, model = "HKY", tstv_ratio = 0, base_frequencies = c(0.3916040, 0.1693122, 0.1620718, 0.2770120)) +
  #sumstat_tajimas_d()+
  #sumstat_four_gamete()

#check_model(model)
#simulate(model)

# one with feat_selection and one without to compare 

# For schistosomes T Krellin et al 2016 
# Generation time = 0.2 years
# mutation rate 8.1×10−9 per basepair per generation 
# rate(theta) is theta = 4Neu 
# top estimate 2000 individuals per host as N (population size) 
# Ne = ( 4 x N breeding males x N breeding females ) / ( N breeding males + N breeding females )


model <- coal_model(sample_size = 76, loci_number = 1, loci_length = 189) 
model <- model + feat_mutation(0.0000015309, model = "HKY", tstv_ratio = 0, base_frequencies = c(0.3916040, 0.1693122, 0.1620718, 0.2770120), gtr_rates = NA, fixed_number = FALSE) + 
         sumstat_tajimas_d() + 
         sumstat_seg_sites() +
         sumstat_four_gamete() +
         sumstat_sfs() + 
         sumstat_nucleotide_div()

sumstats <- simulate(model, seed = 123)
names(sumstats)
sumstats$tajimas_d
sumstats$seg_sites
sumstats$four_gamete
sumstats$pi
sumstats$sfs

# one population of constant size, one genetic locus for 76 haplotypes
#model <- coal_model(sample_size = 76, loci_number = 1, ploidy = 2) +
 # feat_mutation(rate = par_range("theta", 0.1, 5))

#check_model(model)
#result <- simulate(model, seed = 1234)
#result$sfs

# FIND THIS MODEL
#Hudson RR: Generating samples under a Wright-Fisher neutral model
#of neutral variation. Bioinformatics 2002, 18:337-338.

# Coala - coalescent simulation 
# examples https://cran.r-project.org/web/packages/coala/vignettes/coala-examples.html
# vignette("coala-intro", "coala") for intro



model <- coal_model(sample_size = 10, loci_number = 1) +
  feat_mutation(5) +
  sumstat_sfs()
result <- simulate(model)
result$sfs


#neutral distribution of tajimas D, 10 individuals
model <- coal_model(10, 2000) +   # 10 individuals
  feat_mutation(5) +   #5 mutations
  feat_recombination(0) + # no recombination
  sumstat_tajimas_d()
stats <- simulate(model, seed = 15)
plot(density(stats$tajimas_d, na.rm = TRUE), 
     main = "Neutral Distribution of Tajima's D")


# site frequency spectrum - 5 populations and migration
model2a <- coal_model(c(10, 10, 10, 10, 10), 2000) +
  feat_mutation(5) +   # 5 mutations
  feat_recombination(0) +  # no recombination 
  feat_migration(0.5, symmetric = TRUE) +
  sumstat_sfs(population = "all")
stats <- simulate(model2a, seed = 20)
barplot(stats$sfs / sum(stats$sfs), 
        names.arg = seq_along(stats$sfs), 
        col = 3)


# site frequency spectrum - 1 population and no migration
model2a <- coal_model(76, 100) +
  feat_mutation(5) +   # 5 mutations
  feat_recombination(0) +  # no recombination 
  #feat_migration(0, symmetric = TRUE) +
  sumstat_sfs(population = "all")
stats <- simulate(model2a, seed = 20)
barplot(stats$sfs / sum(stats$sfs), 
        names.arg = seq_along(stats$sfs), 
        col = 3, main = "SFS 1 pop - no migration, no recomb, 5 mutations")


# site frequency spectrum but with bottle neck in one population
model2b <- model2a +
  feat_size_change(0.1, population = 2, time = 0.25) +
  feat_size_change(1, population = 2, time = 0.5)
stats <- simulate(model2b, seed = 25)
barplot(stats$sfs / sum(stats$sfs), 
        names.arg = seq_along(stats$sfs), 
        col = 4)



# Plot of nucleotide diversity vs mutation rate 
model3 <- coal_model(10, 50) +
  feat_mutation(par_prior("theta", sample.int(100, 1))) +
  sumstat_nucleotide_div()

stats <- simulate(model3, nsim = 40)
mean_pi <- sapply(stats, function(x) mean(x$pi))
theta <- sapply(stats, function(x) x$pars[["theta"]])
plot(theta, mean_pi, pch = 19, col = "orange")
abline(lm(mean_pi ~ theta), col = "red2", lty = 3)



