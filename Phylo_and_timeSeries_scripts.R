# Scripts for Thuy et al. (2020)
# Miniaturization during a Silurian environmental crisis generated the modern brittle star body plan. Communications Biology

# R code written by David F. Wright

# SET THE WORKING DIRECTORY AND SOURCE FUNCTIONS USED IN ANALYSES
# Must have p.file, t.file, 'Geologic_Timescale.txt' file, and the "functions.R" file in the directory

# Change this to your own directory!
setwd("/Users/Davey/Documents/B_star_phylo_analysis")
source("functions.R")

# INSTALL AND/OR LOAD R PACKAGE DEPENDENCIES

packages <- c("ape")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

# load the library
library(ape)

##############################
# TREE-BASED ANALYSES
##############################

# CALCULATE PROBABILITY OF MORHPOTAXA AS SAMPLED ANCESTORS

# Load output from tip-dated analyses, remove burn-in, and scale the branch lengths to units of time
p.file <- read.table("Tip_dated_tree_analysis.nex.p", header = TRUE)
t.file <- read.nexus("Tip_dated_tree_analysis.nex.t")
trees <- timescaleMrBayes(t.file, p.file); trees <- trees[round(0.25*length(trees) + 2):length(trees)]

# Calculate posterior probabilities that Ophiopetagno_paicei & Muldaster_duplantieri are sampled ancestors
k <- samp.Anc(trees)
petagno_prob <- k[which(names(k) == c("Ophiopetagno_paicei"))]
Muldaster_prob <- k[which(names(k) == c("Muldaster_duplantieri"))]

# Print posterior probabilities of sampled ancestors to screen
print(c(petagno_prob, Muldaster_prob))


# CALCULATE FOSSIL COMPLETNESS THROUGH TIME & Pr[Ancestor-Descendant pairs] per-interval

# load geological timescale and get midpoint times

geol_time<-read.table("Geologic_Timescale.txt",header=TRUE)
System<-cbind(geol_time[,2],geol_time[,3])
row.names(System)<-geol_time[,1]
midpoints <- (System[,1] + System[,2]) / 2

# calculate completeness (Phi), where Phi = (Psi / (Psi + Mu))

foo <- get_Completeness(p.file, param); completeness <- foo$completeness; SD <- foo$SD

# plot completeness through time, default is Paleozoic

Completeness_Through_Time(midpoints, completeness, SD, up_to = 251)

AD_pair <- get_int_AD_pair_probs(p.stat); print(AD_pair)

# Interval with greatest probability of sampling is the Silurian 

AD_pair[which(AD_pair == max(AD_pair))]




##############################
# TIME-SERIES ANALYSES
##############################

# read in data
data <- read.table("Size_data.txt", header = TRUE)

occ_times <- data$Time; Nsamples <- data$Sample_size; M <- data$Plate_area; SD <- data$St_Dev


# PARAMETRIC BOOSTRAP TO GET 95% CONFIDENCE INTERVAL FOR BODY SIZE AT EACH LOCALITY

# parametric boostrap to plot the observed data with 95% confidence interval
sim.values <- sim.Gotland.Series(Nsim = 1000, occ_times, Nsamples, M, SD, bootstrap.M = FALSE, randomize.SD = FALSE)

SE <- get_SE(sim.values$sim.means)

plot(occ_times, M, ylab = expression(Arm ~ plate ~ area ~ (mm)^2), xlab = "Time (Ma)", ylim = c(min(M-SE)-.05 , max(M+SE)+.05), xlim = c(max(occ_times)+1, min(occ_times)-1))

lines(occ_times,M)

points(occ_times,M, cex = 1.2, pch = 19, col = "black"); 

segments(occ_times, M, occ_times, M+2*SE)
segments(occ_times, M, occ_times, M-2*SE);


# MONTE CARLO SIMULATIONS TO GET PROBABILITY OF OBSERVING EXTINCTION PATTERNS BY CHANCE

# 1: tests the magnitude of body size decrease and compares it to a null distribuiton
# 2: tests the directionality of body size decrease for the Mulde Event & Lau Event

# simulate fossil records to use for the two tests:

set.seed(2020) # set seed so figures/histograms can be reproduced exactly
Nsim <- 10000 # 
sim.values <- sim.Gotland.Series(Nsim, occ_times, Nsamples, M, SD, bootstrap.M = TRUE, randomize.SD = FALSE)


# test magnitude of change per extinction event, i.e, whether change from begining to end of each extinction event differs from a null distributon

magnitude.test(sim.values) # Ireviken p = 0.0003, for Mulde and Lau p < 0.00001

# test for the directionality of sustained size decrease across Mulde and Lau Events

# this test is very difficult to reject the null--must have ~large sustained changes
directional.test(sim.values) # Mulde Event p = 0.18 (NS); Lau Event p = 0.05


