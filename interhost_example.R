library(PopGenome)
library(coala)
library(poppr)
require(ape)
library(seqinr)
library(pegas)

# polymorphism here is any substitution or mutation that is present in more than one individual, 
# i.e anything but a singleton.


# Read data into popgenome #
GENOME.class <- readData("/home/milly/Desktop/CD63/fasta", format = "fasta", outgroup = FALSE)

# Available statistics and examples 
show.slots(GENOME.class) 
get.sum.data(GENOME.class)

# Get codon changes - substitution sites, synonymous changes and non synonymous
GENOME.class@region.data@codons
GENOME.class@region.data@synonymous
GENOME.class@region.data@minor.alleles
GENOME.class@region.data@n.singletons

# run nucelotide diveristy and neutrality statistics 
GENOME.class  <- diversity.stats(GENOME.class)
GENOME.class  <- neutrality.stats(GENOME.class)
GENOME.class <- F_ST.stats(GENOME.class)
GENOME.class <-F_ST.stats.2(GENOME.class)

# Detail stats, linkage and sweeps 
GENOME.class <- detail.stats(GENOME.class)
GENOME.class <- linkage.stats(GENOME.class)
GENOME.class <- sweeps.stats(GENOME.class)



# Get nucleotide diveristy and FST statistics 
get.diversity(GENOME.class)[[1]]
GENOME.class@Pi

# Get neutrality and Tajimas D statistic #
get.neutrality(GENOME.class)[[1]]
GENOME.class@Tajima.D
GENOME.class@Fu.Li.F
GENOME.class@Fu.Li.D
GENOME.class@n.segregating.sites
GENOME.class@Fu.F_S
GENOME.class@Strobeck.S

# Ramos and Rozas R_2
GENOME.class@Rozas.R_2

# Linkage 
GENOME.class@Rozas.ZA
GENOME.class@Wall.B
GENOME.class@Wall.Q


# Hudson FST calculations 
GENOME.class@Hudson.G_ST
GENOME.class@Hudson.H_ST
GENOME.class@Hudson.K_ST

# Haplotype counts - n times individuals sequence is in population (n alleles?)
hap_counts<-GENOME.class@region.stats@haplotype.counts

# Recombination test (Hudson.Kaplan.RM) 
GENOME.class <- recomb.stats(GENOME.class)
get.recomb(GENOME.class)[[1]]

#RM:Minimum number of recombination events (Hudson)
GENOME.class@RM

# seg sites = 4 here
GENOME.class@n.segregating.sites

# ape to find base frequencies #
base.freq(cd63_all_data)


# sliding window transform 
GENOME.class.slide<- sliding.window.transform(GENOME.class,width=3,jump=3,type=2,whole.data=TRUE)
GENOME.class.slide<-F_ST.stats(GENOME.class.slide)
GENOME.class.slide <- recomb.stats(GENOME.class.slide)
GENOME.class.slide<-neutrality.stats(GENOME.class.slide)
GENOME.class.slide@Tajima.D
GENOME.class.slide@Pi


# plot sliding window as bar chart for taj d regions...
xaxis <- strsplit(GENOME.class.slide@region.names,split=" ; ")
Pi <- GENOME.class.slide@Pi
taj <- GENOME.class.slide@Tajima.D

#PopGplot(taj) #,colors=TRUE,span=0.1,ylab="Tajima's D",xlab="Codon",
         #ylim=c(min(values,na.rm=TRUE),max(values,na.rm=TRUE)))

# work on this as 0 means somethign for taj d 
tajd<-data.frame(GENOME.class.slide@Tajima.D)
tajd[is.na(tajd)] <- 0
PopGplot(tajd)


# Define populations #
GENOME.class_interhost <- set.populations(GENOME.class, list(c("12b", "12c", "12d", "12e", "12f", "12g"), c("13a", "13c", "13f"), c("14a", "14b","14c","14d"), c("15a","15b","15c"),c("16c","16d","16e"),c("20a","20b","20d"), c("A01","A02","A03","A04","A05", "A06","A07","A08","A09","A10", "A11", "A12", "A13"), c( "B01","B02","B03", "B04", "B05", "B06", "B07","B08","B09", "B10", "B12"), c("CA","CC","CD","CE","CF","CG","CI","CJ","CL","CM"), c("DB",  "DC" , "DD" ,"DF"  ,"DG" , "DH",  "DI", "DK",  "DL",  "DM"), c("EA","EB" , "EC",  "EE" , "EF",  "EG",  "EI",  "EJ", "EL",  "EM")))
GENOME.class_interhost  <- neutrality.stats(GENOME.class_interhost)
GENOME.class_interhost@Tajima.D
GENOME.class_interhost@n.biallelic.sites
GENOME.class_interhost <-diversity.stats(GENOME.class_interhost)
GENOME.class_interhost <-F_ST.stats(GENOME.class_interhost)
GENOME.class_interhost@Pi

## need to be more than 4 individuals #
GENOME.class_interhost@Fu.Li.F
GENOME.class_interhost@Fu.Li.D
GENOME.class_interhost@Fu.F_S

## pairwise fst comparisons ##
GENOME.class_interhost@nuc.F_ST.pairwise
GENOME.class_interhost@nuc.F_ST.vs.all


get.sum.data(GENOME.class_interhost)

GENOME.class_interhost@Fu.F_S

## hapltype diversity within and between populations ##
GENOME.class_interhost@hap.diversity.between 
GENOME.class_interhost@hap.diversity.within

##nucleotide diveristy between and within ##
GENOME.class_interhost@nuc.diversity.between
GENOME.class_interhost@nuc.diversity.within


