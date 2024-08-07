#####################
# PlastidHaplotypes.R
# Jessie Pelosi
# Last Modified July 13, 2022
#
#####################

########## Generate Haplotype Networks ########## 
# Adapted from https://cran.r-project.org/web/packages/pegas/vignettes/PlotHaploNet.pdf
library(ape)
library(pegas)
library(haplotypes)

# Path for individual-level subset: PopulationGenetics/RADOrgMiner/80subset/
# Path for colony-level data: PopulationGenetics/RADOrgMiner/Colonies/

popmap <- read.delim("../PopulationGenetics/RADOrgMiner/Colonies/PopMap.txt", header= F)
plastid <- read.FASTA("../PopulationGenetics/RADOrgMiner/Colonies/concat_loci.fa")
plastid.dna<- as.dna(plastid)
plastid_nogaps <- haplotypes::remove.gaps(x= plastid.dna, entire.col=T)
plastid_nogaps.bin <- as.DNAbin(plastid_nogaps)
plastid_haplo <- pegas::haplotype(plastid_nogaps.bin)
length(seg.sites(plastid_nogaps.bin))
summary(plastid_haplo)
plastid_net <- haploNet(plastid_haplo)
plot(plastid_net, size = attr(plastid_net, "freq"), fast = FALSE)
re <- replot()

# Run on populations 
# For colony-level data, popmap$V1, for subset data used popmap$V2
R <- haploFreq(plastid_nogaps.bin, fac=popmap$V1, haplo=plastid_haplo)
size <- summary(plastid_haplo)
dist <- dist.dna(plastid_haplo, "N")
nt <- rmst(dist)
nt.labs <- attr(nt, "labels")
size <- size[nt.labs]
plot(nt, size = size, show.mutation = T)
R <- R[nt.labs, ]

plot(nt, size = size, pie =R, legend = c(-40,70))
re <- replot()

# Run on colonies
# For colony-level data, popmap$V2, for subset data used popmap$V3
C <- haploFreq(plastid_nogaps.bin, fac=popmap$V2, haplo=plastid_haplo)
size <- summary(plastid_haplo)
dist <- dist.dna(plastid_haplo, "N")
nt <- rmst(dist)
nt.labs <- attr(nt, "labels")
size <- size[nt.labs]
plot(nt, size = size, show.mutation = T)
C <- C[nt.labs, ]

plot(nt, size = size, pie =C, legend = c(-50,15))
re <- replot()

nuc.div(plastid, pairwise.deletion = F, variance = T) #On full set = 0.00049
nuc.div(plastid_haplo, pairwise.deletion = F) #On haplotypes = 0.00049
nuc.div(plastid, pairwise.deletion = T) #Remove gaps = 0.00256168
nuc.div(plastid_haplo, pairwise.deletion = T) #On haplotypes, remove gaps = 0.00049

########## Calculate Statistics ########## 
# Adapted from https://popgen.nescent.org/PopDiffSequenceData.html
# And https://grunwaldlab.github.io/Population_Genetics_in_R/Population_Strata.html
library(apex)
library(adegenet)
library(poppr)
library(mmod)

# Read in each locus individually
files <- dir("../PopulationGenetics/RADOrgMiner/80subset/", pattern = "^loc.*.fa", full = T)
loci <- read.multiFASTA(files)

# Replace periods and dashes in locus names (otherwise this throws an error in next command)
getLocusNames(loci)
setLocusNames(loci) <- gsub(".fa", "",getLocusNames(loci))
setLocusNames(loci) <- gsub("-","", getLocusNames(loci))
getLocusNames(loci)

# Concatenate loci and convert to genid object
plastid.gid <- multidna2genind(loci, mlst = TRUE)

# Add strata info
population_info <- read.delim("../PopulationGenetics/RADOrgMiner/80subset/PopMapNum.txt", header = T)
strata <- data.frame(populations=population_info$Pop, colonies=population_info$Colony)
strata(plastid.gid) <- strata
setPop(plastid.gid) <- ~colonies

# Quick basic stats for colonies 
diff_stats(plastid.gid) 

# Hs        Ht   Gst_est Gprime_st     D_het    D_mean 
# 0.2033661 0.5551430 0.6336691 0.8015505 0.4509744 0.3771379 

Phi_st_Meirmans(plastid.gid) #global = 0.1636

plastid.gid_noclone <- clonecorrect(plastid.gid, strata = ~colonies) 
diff_stats(plastid.gid_noclone) # Values are same as above 
Phi_st_Meirmans(plastid.gid_noclone) # Values are same as above 


# Set populations as populations, rather than colonies 
setPop(plastid.gid) <- ~populations

# Quick basic stats for populations
plastid.gcl <- as.genclone(plastid.gid)
divstats <- poppr(plastid.gcl)
divstats
diff_stats(plastid.gid) 

#       Hs        Ht       Gst_est   Gprime_st   D_het    D_mean 
# 0.46444344 0.54653658 0.15020614 0.30008344 0.16605948 0.09314214  

Phi_st_Meirmans(plastid.gid) #global = 0.0603

###############
# Merge colonies and re-run 

#####################
# PlastidHaplotypes.R
# Jessie Pelosi
# Last Modified July 13, 2022
#
#####################

########## Generate Haplotype Networks ########## 
# Adapted from https://cran.r-project.org/web/packages/pegas/vignettes/PlotHaploNet.pdf
library(ape)
library(pegas)
library(haplotypes)

popmap <- read.delim("../PopulationGenetics/RADOrgMiner/Colonies/PopMap.txt", header= F)
plastid <- read.FASTA("../PopulationGenetics/RADOrgMiner/Colonies/concat_loci.fa")
plastid.dna<- as.dna(plastid)
plastid_nogaps <- haplotypes::remove.gaps(x= plastid.dna, entire.col=T)
plastid_nogaps.bin <- as.DNAbin(plastid_nogaps)
plastid_haplo <- pegas::haplotype(plastid_nogaps.bin)
length(seg.sites(plastid_nogaps.bin))
summary(plastid_haplo)
plastid_net <- haploNet(plastid_haplo)
plot(plastid_net, size = attr(plastid_net, "freq"), fast = FALSE)
re <- replot()

#Run on populations 
R <- haploFreq(plastid_nogaps.bin, fac=popmap$V2, haplo=plastid_haplo)
size <- summary(plastid_haplo)
dist <- dist.dna(plastid_haplo, "N")
nt <- rmst(dist)
nt.labs <- attr(nt, "labels")
size <- size[nt.labs]
plot(nt, size = size, show.mutation = T)
R <- R[nt.labs, ]

plot(nt, size = size, pie =R, legend = c(-40,70))
re <- replot()

# Run on colonies
C <- haploFreq(plastid_nogaps.bin, fac=popmap$V3, haplo=plastid_haplo)
size <- summary(plastid_haplo)
dist <- dist.dna(plastid_haplo, "N")
nt <- rmst(dist)
nt.labs <- attr(nt, "labels")
size <- size[nt.labs]
plot(nt, size = size, show.mutation = T)
R <- R[nt.labs, ]

plot(nt, size = size, pie =R, legend = c(-55,15))
re <- replot()

nuc.div(plastid, pairwise.deletion = F, variance = T) #On full set = 0.00049
nuc.div(plastid_haplo, pairwise.deletion = F) #On haplotypes = 0.00049
nuc.div(plastid, pairwise.deletion = T) #Remove gaps = 0.00256168
nuc.div(plastid_haplo, pairwise.deletion = T) #On haplotypes, remove gaps = 0.00049

########## Calculate Statistics ########## 
# Adapted from https://popgen.nescent.org/PopDiffSequenceData.html
# And https://grunwaldlab.github.io/Population_Genetics_in_R/Population_Strata.html
library(apex)
library(adegenet)
library(poppr)
library(mmod)

# Read in each locus individually
files <- dir("../PopulationGenetics/RADOrgMiner/Colonies/", pattern = "^loc.*.fa", full = T)
loci <- read.multiFASTA(files)

# Replace periods and dashes in locus names (otherwise this throws an error in next command)
getLocusNames(loci)
setLocusNames(loci) <- gsub(".fa", "",getLocusNames(loci))
setLocusNames(loci) <- gsub("-","", getLocusNames(loci))
getLocusNames(loci)

# Concatenate loci and convert to genid object
plastid.gid <- multidna2genind(loci, mlst = TRUE)

# Add strata info
population_info <- read.delim("../PopulationGenetics/RADOrgMiner/Colonies/PopMap.txt", header = F)
strata <- data.frame(populations=population_info$V1, colonies=population_info$V2)
strata(plastid.gid) <- strata
setPop(plastid.gid) <- ~colonies

# Quick basic stats for colonies 
diff_stats(plastid.gid, phi_st = T) 

# Hs        Ht   Gst_est Gprime_st     D_het    D_mean 
# 0.38398734 0.46889740 0.18108452 0.31372497 0.14932469 0.08803831 

Phi_st_Meirmans(plastid.gid) #global = -0.045 

plastid.gid_noclone <- clonecorrect(plastid.gid, strata = ~colonies) 
diff_stats(plastid.gid_noclone) # Values are same as above 
Phi_st_Meirmans(plastid.gid_noclone) # Values are same as above 


# Set populations as populations, rather than colonies 
setPop(plastid.gid) <- ~populations

# Quick basic stats for populations
plastid.gcl <- as.genclone(plastid.gid)
divstats <- poppr(plastid.gcl)
divstats
diff_stats(plastid.gid) 

#       Hs        Ht       Gst_est   Gprime_st   D_het    D_mean 
# 0.46444344 0.54653658 0.15020614 0.30008344 0.16605948 0.09314214  

Phi_st_Meirmans(plastid.gid) #global = 0.0603