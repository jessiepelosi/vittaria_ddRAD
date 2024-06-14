#####################
# NuclearddRADAnalyses.R
# Jessie Pelosi
# Last Modified June 2024
#####################

library(adegenet)
library(pegas)
library(hierfstat)
library(dplyr)
library(vcfR)
library(ggplot2)
library(poppr)
library(HWxtest)
library(nlme)
library(dartR)
library(RColorBrewer)

##### -------------------------- Read in data ----------------------------#####

# Read in population assignment/information for the nuclear subsetted data 
population_info <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_PopCol.txt", header = T)

# Read in data for up to 50% missing data
sub50 <- read.vcfR("../PopulationGenetics/ddRADSeq/ipyrad/Individual-Level/Subset_95Clust/Subset_Clust95_filt0.5_bi.recode.vcf", 
                   convertNA = T)

sub50.gid <- vcfR2genind(sub50, strata = population_info, return.alleles = T)
sub50.gid
pop(sub50.gid) <- population_info$Population
head(strata(sub50.gid))

# Write basyescan file 

data_fstat <- genind2hierfstat(sub50.gid)
write.bayescan(dat=data_fstat, diploid=T, fn="sub50.bsc")

##### ------------------------- Allelic Richness ---------------------------#####

sum <- summary(sub50.gid)
pop(sub50.gid) <- population_info$Population
popPA <- as.data.frame(private_alleles(sub50.gid) %>%  apply(MARGIN =1, FUN = sum))

pop(sub50.gid) <- population_info$Colony
colPA <- as.data.frame(private_alleles(sub50.gid) %>%  apply(MARGIN =1, FUN = sum))

# Rareified Allelic Richness
pop(sub50.gid) <- population_info$Population
sub50.loci <- as.loci(sub50.gid)
PopAr <- pegas::allelicrichness(sub50.loci, method = "extrapolation", min.n = NULL, pop = sub50.loci$population)
PopAr_melt <- reshape2::melt(PopAr)
plot(x = PopAr_melt$Var1, y = PopAr_melt$value)
PopArsummary <- PopAr_melt %>% 
  group_by(Var1) %>% 
  summarize(mean = mean(value))

min(PopArsummary$mean); mean(PopArsummary$mean); max(PopArsummary$mean)
# Min: 1.198376; Mean: 1.337151; Max: 1.487612

pop(sub50.gid) <- population_info$Colony
sub50.loci <- as.loci(sub50.gid)
ColonyAr <- pegas::allelicrichness(sub50.loci, method = "extrapolation", min.n = NULL, pop = sub50.loci$population)
ColonyAr_melt <- reshape2::melt(ColonyAr)
plot(x = ColonyAr_melt$Var1, y = ColonyAr_melt$value)
ColonyARsummary <- ColonyAr_melt %>% 
  group_by(Var1) %>% 
  summarize(mean = mean(value))

min(ColonyARsummary$mean); mean(ColonyARsummary$mean); max(ColonyARsummary$mean)
# Min: 1.004588; Mean: 1.115536; Max:1.28524


##### ------------------------ Basic popgen stats --------------------------#####

pop(sub50.gid) <- population_info$Population
basic.stats(sub50.gid)
# Population-level (population = populations)
#   Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.0747  0.0736  0.0920  0.0184  0.0936  0.0200  0.1998  0.2134 -0.0148  0.0216

pop(sub50.gid) <- population_info$Colony
basic.stats(sub50.gid)
# Colony-level (population = colony)
#   Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.0721  0.0680  0.0857  0.0177  0.0863  0.0183  0.2065  0.2116 -0.0598  0.0196 

# Expected heterozygosity at population-level
Hs(data = sub50.gid)  
#     DB         PR         CC         PN         RP         TG         NF   
# 0.07125675 0.06980915 0.07757789 0.05013127 0.09068627 0.06309888 0.09530340
#     RC         PM         BP        AA         MM         SH 
# 0.07786984 0.06884393 0.09005882 0.06824656 0.10660981 0.06738062 

# Observed heterozygosity at population-level 
Ho_pop <- as.data.frame(Ho(data =sub50.gid)) 
#     DB         PR         CC         PN         RP         TG         NF        
#0.07562692 0.07683756 0.07090205 0.06510990 0.09798808 0.07060494 0.09995683 
#    RC         PM         BP         AA         MM         SH 
#0.07327378 0.07890124 0.07476669 0.07552502 0.07183908 0.06999163 

# Observed heterozygosity at colony-level
pop(sub50.gid) <- population_info$Colony
Ho_col <- as.data.frame(Ho(data = sub50.gid))  
min(Ho_col$`Ho(data = sub50.gid)`)
mean(Ho_col$`Ho(data = sub50.gid)`)
max(Ho_col$`Ho(data = sub50.gid)`)

# Expected heterozygosity at colony-level
Hs(data = sub50.gid)
Hs_col <- as.data.frame(Hs(data = sub50.gid))  
min(Hs_col$`Hs(data = sub50.gid)`, na.rm = T)
mean(Hs_col$`Hs(data = sub50.gid)`, na.rm = T)
max(Hs_col$`Hs(data = sub50.gid)`, na.rm = T)

#####-------------------------- F Statistics --------------------------######

WC <- wc(sub50.gid) 
# Populations (population = population)           Colonies (population = colony)
#   FST 0.1962578                                  FST  0.2284858
#   FIS -0.002418957                               FIS -0.0530314

pfst <- pairwise.WCfst(sub50.gid)
btpft <- boot.ppfst(sub50.gid)

# Get proportion of loci that have negative Fis values for each population
pop(sub50.gid) <- population_info$Population
obj <- seppop(sub50.gid)
Fis <- basic.stats(obj$DB); tmp <- as.data.frame(Fis$perloc$Fis)
filter(tmp, Fis$perloc$Fis < 0) %>% dplyr::summarise(n=n()); filter(tmp, Fis$perloc$Fis > 0) %>% dplyr::summarise(n=n())

boot.ppfis(sub50.gid)

pop(sub50.gid) <- population_info$Colony
obj <- seppop(sub50.gid)

basic.stats(obj$TG3)

##### ----------------------------- AMOVA ------------------------------- #####

# Convert vcf to genlight object, poppr.amova doesn't handle genind files well
# as poppr treats zerso differently to handle mixed or ambiguous ploidy 
sub50.gl <- vcfR2genlight(sub50)
pop(sub50.gl) <- population_info$Population
ploidy(sub50.gl) <- 2
strata(sub50.gl) <- population_info

# With ade4 method 
#sub50.miss <- missingno(sub50.gid, type = "mean")
amova.res <- poppr.amova(sub50.gl, ~Population/Colony, method = "ade4") 
amova.res
amova.test <- randtest(amova.res, nrepet = 999) 
plot(amova.test)
amova.test 

##### ---------------------------- Linkage ------------------------------ #####

# Total rd for all 82 samples 
ia(sub50.gid, sample = 999) 
#     Ia        p.Ia       rbarD        p.rD                                                                     
#23.21557223  0.00100000  0.01744396  0.00100000 

# Get rd for each population
rd_test <- sapply(seppop(sub50.gid), 
               function(ls) ia(ls, sample = 999))

# Get rd for each colony 
pop(sub50.gid) <- population_info$Colony
rd_test <- sapply(seppop(sub50.gid), 
                  function(ls) ia(ls, sample = 999))

# Compare per-site pi values 
total_pi <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/out.sites.pi")
mean(total_pi$PI) #average total pi=0.0869909 
AA_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/AA.sites.pi"); AA_pi <- mean(AA_pi_df$PI, na.rm = T) 
BP_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/BP.sites.pi"); BP_pi <- mean(BP_pi_df$PI, na.rm = T) 
CC_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/CC.sites.pi"); CC_pi <- mean(CC_pi_df$PI, na.rm = T)
DB_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/DB.sites.pi"); DB_pi <- mean(DB_pi_df$PI, na.rm = T)
MM_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/MM.sites.pi"); MM_pi <- mean(MM_pi_df$PI, na.rm = T)
NF_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/NF.sites.pi"); NF_pi <- mean(NF_pi_df$PI, na.rm = T)
PM_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/PM.sites.pi"); PM_pi <- mean(PM_pi_df$PI, na.rm = T)
PN_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/PN.sites.pi"); PN_pi <- mean(PN_pi_df$PI, na.rm = T)
PR_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/PR.sites.pi"); PR_pi <- mean(PR_pi_df$PI, na.rm = T)
RC_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/RC.sites.pi"); RC_pi <- mean(RC_pi_df$PI, na.rm = T)
RP_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/RP.sites.pi"); RP_pi <- mean(RP_pi_df$PI, na.rm = T)
SH_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/SH.sites.pi"); SH_pi <- mean(SH_pi_df$PI, na.rm = T)
TG_pi_df <- read.delim("../PopulationGenetics/ddRADSeq/ipyrad/Subset_95Clust/TG.sites.pi"); TG_pi <- mean(TG_pi_df$PI, na.rm = T)
pop_names <- c("AA", "BP", "CC", "DB", "MM", "NF", "PM", "PN", "PR", "RC", "RP", "SH", "TG")
pi <- c(AA_pi, BP_pi, CC_pi, DB_pi, MM_pi, NF_pi, PM_pi, PN_pi, PR_pi, RC_pi, RP_pi, SH_pi, TG_pi)
pi_df <- data.frame(pop_names, pi)
ggplot(data = pi_df, mapping = aes(x=pop_names, y=pi)) + geom_col()

##### --------------------------- Clustering --------------------------------#####
tre <- nj(pfst)
plot(tre, type = "unrooted") 

grp <- find.clusters(sub50.gid, max.n.clust = 50)
dapc <- dapc(sub50.gid, grp$grp)
scatter(dapc)
grp
compoplot(dapc)

# Re-draw structure-like plot 
# https://luisdva.github.io/rstats/dapc-plot/
library(ggh4x)
postprobs <- as.data.frame(round(dapc$posterior, 4))
clusters <- tibble::rownames_to_column(postprobs, var = "ind") %>% 
  mutate(site = population_info$Population, colony = population_info$Colony)
clusters_long <- reshape2::melt(clusters)
ggplot(data = clusters_long, mapping = aes(x = ind, y = value, fill = factor(variable))) + 
  geom_col(color = "gray", size=0.01) + 
  facet_nested(~site + colony, switch = "x", nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free") +
  theme_minimal() +
  labs(x = "Samples", y = "Membership Probability") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  theme(panel.spacing.x = unit(0.18, "lines"), 
        axis.text.x = element_blank(), panel.grid = element_blank())

ggsave("DPAC_structure_plot.pdf", height = 5, width = 30)

snapclust.choose.k(13, sub50.gid) # k=2 is best fit with snapclust 
res <- snapclust(sub50.gid, k =2) 
compoplot(res)

######--------------------------- Assess clonality ---------------------------#####

## Explore data clustering  

tree <- aboot(sub50.gl, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
tree
nb.cols <- 13
cols <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
plot.phylo(tree, cex = 0.6, font = 2, adj = 0, tip.color =  cols[pop(sub50.gl)])
nodelabels(tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)

library(igraph)
sub50.dist <- bitwise.dist(sub50.gl)
sub50.msn <- poppr.msn(sub50.gl, sub50.dist, showplot = FALSE, include.ties = T)
node.size <- rep(2, times = nInd(sub50.gl))
names(node.size) <- indNames(sub50.gl)
vertex.attributes(sub50.msn$graph)$size <- node.size
set.seed(9)
plot_poppr_msn(sub50.gl, sub50.msn , palette = cols, gadj = 70)

sub50.pca <- glPca(sub50.gl, nf = 3)
barplot(100*sub50.pca$eig/sum(sub50.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

sub50.pca.scores <- as.data.frame(sub50.pca$scores)
sub50.pca.scores$pop <- pop(sub50.gl)
p <- ggplot(sub50.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p

## Generate multi-locus lineages 

# Calculate pairwise distances between samples  
# provesti.dist is the same as diss.dist, just as % not raw allelic differences 
distance <- provesti.dist(sub50.gid)
distance
write.table(as.matrix(distance), file = "provesti_dist.txt")
PopDiff <- reshape2::melt(read.delim("ProvestiPop.txt", sep = ''))
ggplot(data = PopDiff, mapping = aes(x = value, fill = variable)) + geom_density(alpha = 0.5, position = "identity") +
  theme_classic() + xlab("Genetic Distance")
t.test(value ~ variable, data = PopDiff, alternative = "less")
ggsave("provesti_population_distance.png", dpi = 300, height = 5, width =7)


ColDiff <- reshape2::melt(read.delim("ProvestiCol.txt", na.strings = "0"))
ggplot(data = ColDiff, mapping = aes(x = value, fill = variable)) + geom_density(alpha = 0.5, position = "identity") +
  theme_classic() + xlab("Genetic Distance")
t.test(value ~ variable, data = ColDiff, alternative = "less")
ggsave("provesti_colony_distance.png", dpi = 300, height = 5, width =7)

sub50_filtered <- filter_stats(sub50.gc, distance = distance, plot = TRUE)

# Alter threshold values to explore clonality and diversity in colonies 

# Threshold = 0.004 #
mlg.filter(sub50.gc, distance = distance, algorithm = "farthest_neighbor") <- 0.004
sub50.gc
thresh0.004 <- mlg.table(sub50.gc)

sub50.cc<- clonecorrect(sub50.gc, strata = ~Population/Colony)
basic.stats(sub50.cc)
#   Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.0766  0.0759  0.0947  0.0188  0.0963  0.0204  0.1986  0.2121 -0.0093  0.0221 
sub50.cc.miss <- missingno(sub50.cc, type = "mean")
poppr.amova(sub50.cc.miss, ~Population/Colony, algorithm = "ade4")

# Threshold = 0.01 #
mlg.filter(sub50.gc, distance = distance, algorithm = "farthest_neighbor") <- 0.01
sub50.gc
thresh0.01 <- mlg.table(sub50.gc)

sub50.cc <- clonecorrect(sub50.gc, strata = ~Population/Colony)
basic.stats(sub50.cc)
#   Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
#0.0765  0.0759  0.0949  0.0190  0.0965  0.0206  0.1999  0.2135 -0.0084  0.0223 
sub50.cc.miss <- missingno(sub50.cc, type = "ignore")
poppr.amova(sub50.cc.miss, ~Population/Colony)

# Threshold = 0.02 #
mlg.filter(sub50.gc, distance = distance, algorithm = "farthest_neighbor") <- 0.02
sub50.gc
thresh0.02 <- mlg.table(sub50.gc)

sub50.cc <- clonecorrect(sub50.gc, strata = ~Population/Colony)
basic.stats(sub50.cc)
#   Ho     Hs     Ht    Dst    Htp   Dstp    Fst   Fstp    Fis   Dest 
#0.0757 0.0762 0.0946 0.0184 0.0963 0.0201 0.1948 0.2083 0.0072 0.0217 

sub50.cc.miss <- missingno(sub50.cc, type = "mean")
poppr.amova(sub50.cc.miss, ~Population/Colony)

######---------------------- Within-colony diversity ---------------------- ######

PN2 <- read.vcfR("within_colonies_sub50/PN2.recode.vcf", convertNA = T)
PN2.gid <- vcfR2genind(PN2); PN2.loci <- genind2loci(PN2.gid)
distgenDIFF <- dist.gene(PN2.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
distgenDISS <- diss.dist(PN2.gid, percent = TRUE, mat = FALSE)
pop(PN2.gid) <- c("PN2", "PN2", "PN2"); basic.stats(PN2.gid)


CC9 <- read.vcfR("within_colonies_sub50/CC9.recode.vcf", convertNA = T)
CC9.gid <- vcfR2genind(CC9); CC9.loci <- genind2loci(CC9.gid)
distgenDIFF <- dist.gene(PN2.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
distgenDISS <- diss.dist(PN2.gid, percent = TRUE, mat = FALSE)
pop(CC9.gid) <- c("CC9", "CC9", "CC9"); basic.stats(CC9.gid)

RC1 <- read.vcfR("within_colonies_sub50/RC1.recode.vcf", convertNA = T)
RC1.gid <- vcfR2genind(RC1); RC1.loci <- genind2loci(RC1.gid)
distgenDIFF <- dist.gene(RC1.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
distgenDISS <- diss.dist(RC1.gid, percent = TRUE, mat = FALSE)
pop(RC1.gid) <- c("RC1", "RC1", "RC1", "RC1"); basic.stats(RC1.gid)

PM2 <- read.vcfR("within_colonies_sub50/PM2.recode.vcf", convertNA = T)
PM2.gid <- vcfR2genind(PM2); PM2.loci <- genind2loci(PM2.gid)
distgenDIFF <- dist.gene(PM2.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
distgenDISS <- diss.dist(PM2.gid, percent = TRUE, mat = FALSE)
pop(PM2.gid) <- c("PM2", "PM2", "PM2", "PM2"); basic.stats(PM2.gid)

AA3 <- read.vcfR("within_colonies_sub50/AA3.recode.vcf", convertNA = T)
AA3.gid <- vcfR2genind(AA3); AA3.loci <- genind2loci(AA3.gid)
distgenDIFF <- dist.gene(AA3.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
distgenDISS <- diss.dist(AA3.gid, percent = TRUE, mat = FALSE)
pop(AA3.gid) <- c("AA3", "AA3", "AA3"); basic.stats(AA3.gid)

SH7 <- read.vcfR("within_colonies_sub50/SH7.recode.vcf", convertNA = T)
SH7.gid <- vcfR2genind(SH7); SH7.loci <- genind2loci(SH7.gid)
distgenDIFF <- dist.gene(SH7.loci, method="pairwise", pairwise.deletion = FALSE, variance = FALSE)
distgenDISS <- diss.dist(SH7.gid, percent = TRUE, mat = FALSE)
pop(SH7.gid) <- c("SH7", "SH7", "SH7"); basic.stats(SH7.gid)

win_col <- read.csv("within_colony_mantel.csv")
ggplot(data = win_col, mapping = aes(x = Physical.Distance, y = Ratio.Genetic.Distance)) +
  geom_point() + theme_classic() + geom_smooth(method = "lm") + facet_wrap(~Colony)

ggsave("within_colony_mantel_facet.png", dpi = 300, height = 7, width = 10)

lin <- nlme::lme(Ratio.Genetic.Distance ~ Physical.Distance, random = ~1|Colony, data = win_col, na.action = na.omit)
summary(lin)

# Linear mixed-effects model fit by REML
#  Data: win_col 
#      AIC       BIC     logLik
#  -88.08222 -84.09929 48.04111
#
# Random effects:
#  Formula: ~1 | Colony
#  (Intercept)   Residual
#  StdDev: 1.055537e-06 0.01598041
#  Fixed effects:  Ratio.Genetic.Distance ~ Physical.Distance 
#                      Value      Std.Error DF  t-value p-value
#  (Intercept)       0.04305030 0.005328147 15 8.079788  0.0000
#  Physical.Distance 0.00016797 0.000136550 15 1.230133  0.2376
#  Correlation: 
#  (Intr)
#  Physical.Distance -0.769
#
# Standardized Within-Group Residuals:
#  Min          Q1         Med          Q3         Max 
#  -1.16553705 -0.90397472 -0.05046153  0.53414803  2.25872792 
#  Number of Observations: 22
#  Number of Groups: 6 


### ----------- Bayescan ------------- ### 

bs <- read.table("sub50_fst.txt")

filter(bs, qval < 0.05)



## Plot dadi model vs empirical fs

table <- read.delim("threeEpoch_values.txt")

table$y2 <- as.numeric(table$y2)

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

ggplot(data = table, mapping = aes(x = x, y = y2, color = FSTYPE)) + geom_point() + geom_line() + 
  xlim(0, 25) + scale_y_continuous(trans = log10_trans(), breaks=base_breaks()) + theme_classic() + 
  xlab("Frequency Class of Alleles") + ylab("Number of Alleles") + scale_color_manual(values = c("fs" = "dodgerblue3", "model_fs"= "firebrick2"))

ggsave("sfs-model-comp.png", dpi=300, height = 5, width = 9)
