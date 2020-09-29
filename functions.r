# Library of phylo-related functions (dependent on 'ape') that might be useful here and some that just streamline viewing the output for very specific analysis/plots

# Functions to simulate brittle star size evolution and Monte Carlo analyses to obtain probabilities of extinction patterns begin on line 363

# by David F. Wright (We are 138!)

# Function to get posterior probabilities of a taxon being a sampled ancestor
# Note the the "trees" argument must be a post-burn-in sample of the posterior distribution

samp.Anc <- function(trees){
  
  freq <- vector(mode = "numeric", length = Ntip(trees[[1]]))
  names(freq) <- trees[[1]]$tip.label
  
  for (i in 1:length(trees)){
    
    tip.br <- setNames(trees[[i]]$edge.length[sapply(1:length(trees[[i]]$tip.label),function(x,y)   which(y == x),y = trees[[i]]$edge[,2])],trees[[i]]$tip.label)
    
    for (j in 1:length(tip.br)){
      
      if (tip.br[j] == 0){
        freq[j] <- freq[j] + 1
        
      }
      
    }
    
  }
  
  anc.prop <- (freq / length(trees))
  return(anc.prop)
  
}

# FUNCTION TO CHANGE BRANCH LENGTHS BELOW A THRESHOLD TO A CONSTANT  
# Can take either a list of trees or a single tree. Can also be used to scale branches to unit length when NULL (e.g., a parsimony tree estimated in PHANGORN, etc.)

modify.edges <- function(trees, blens = 1, threshold = 1e-01){
  # Check whether object 'trees' is a single tree or list of trees
  # For multi-phylo objects
  
  if (class(trees) == "multiPhylo"){
    
    for(i in 1:length(trees)){  # loop over all trees a list of trees
      
      tree <- trees[[i]]
      
      # check to make sure there are branch lengths; otherwise add small value to get rid of edge.length == NULL warnings()
      if (is.null(tree$edge.length)){
        N.edges <- (Ntip(tree) - 3) + Ntip(tree) # == number of internal branches + number of tips
        
        for (j in 1:N.edges){
          
          tree$edge.length[j] <- 1e-02
          
        }
        
      }
      
      # Main loop to add edge.lengths below threshold for each jth tree
      for (j in 1:length(tree$edge.length)){
        
        edge <- tree$edge.length[j]
        
        if (edge <= threshold){
          tree$edge.length[j] <- blens
          
        }
        
      }
      
      trees[[i]] <- tree
      
    }
    
    return(trees) # returns a list of trees with modified edges
    
  }
  
  # For a single tree
  if (class(trees) == "phylo"){
    tree <- trees
    
    # check to make sure there are branch lengths; otherwise add small value to get rid of edge.length == NULL warnings()
    if (is.null(tree$edge.length)){
      N.edges <- (Ntip(tree) - 3) + Ntip(tree) # == number of internal branches + number of tips
      
      for (j in 1:N.edges){
        
        tree$edge.length[j] <- 1e-02
        
      }
      
    }
    
    # Main loop to add edge.lengths below threshold 
    for (k in 1:length(tree$edge.length)){
      
      edge <- tree$edge.length[k]
      
      if (edge <= threshold){
        tree$edge.length[k] <- blens
        
      }
      
    }
    
    return(tree) # returns a single tree with modified edges
    
  }
  
}


# Convert MrBayes estimates of diversification dynamics & fossilization to speciation, extinction, and sampling rates

calcMu <- function(netDiv, turnover){
  
  mu <- (netDiv * turnover) / (1 - turnover)
  return(mu)
  
}

calcLambda <- function(netDiv, mu){
  
  lambda <- netDiv + mu
  return(lambda)
  
}

calcPsi <- function(mu, relSamp){
  
  Psi <- (mu * relSamp) / (1 - relSamp)
  return(Psi)
  
}


# Function to add length to Zero Length Branches
removeZLB <- function(trees, blens = 1){
  
  for(i in 1:length(trees)){
    
    tree <- trees[[i]]
    
    for (j in 1:length(tree$edge.length)){
      
      edge <- tree$edge.length[j]
      
      if (edge <= 10e-6){
        tree$edge.length[j] <- blens
        
      }
      
    }
    
    trees[[i]] <- tree
    
  }
  
  return(trees)
  
}


# Sample MrBayes trees from the posterior and time scale them in units of relative time (must use anchorTree function to scale to absolute time)
sample_timescaleMrBayes <- function(t.file, p.file, subsample = 100, burnin = 0.35){ 
 
  n <- length(t.file)
  clockrate <- p.file$clockrate
  samples <- sample(seq(from = round(n*burnin), to=n), size = subsample)
  scaled.trees <- list()
  
     for (i in 1:length(samples)) {
  
        tree <- samples[i]
        foo <- t.file[[tree]]
        foo$edge.length <-  foo$edge.length / clockrate[tree] 
        scaled.trees[[i]] <- foo
  
     }
 
  return(scaled.trees)

}



# Scales all trees in the t.file to relative time. Use anchorTree to scale the tree to absolute time after using this function.
timescaleMrBayes <- function(t.file, p.file){
  
  for (i in 1:length(t.file)){
  
    clock.rate <- p.file$clockrate[i]
    t.file[[i]]$edge.length <- t.file[[i]]$edge.length / clock.rate
  
  }
 
  return(t.file)

}

# Find a tree in the posterior distribution t.file. May be useful to get associated info from p.file for a given tree (e.g., MCCT). May be slow if #trees is large.
findTree <- function(tree, t.file){
  
  clades <- vector(length = length(t.file))
 
   for (i in 1:length(t.file)){
    
     clades[i] <- dist.topo(t.file[[i]], tree)
  
   }
  
  return(which(clades == 0)) 

} 


# get distribution for the age of the age of the clade

getCladeAge <- function(trees){
  
  age <- c()
  
  for (i in 1:length(trees)){
    
    age[i] <- max(dateNodes(trees[[i]]))
    
  }
  
  return(age)
  
}


  # Use anchor taxon to scale the tree to abosulte time. This only makes sense if trees are already
  # scaled to relative time
  anchorTree <- function(trees, anchorTime){
    # for a list of trees
    if (class(trees) == "multiPhylo"){
      for (i in 1:length(trees)){
        oldRootAge <- max(dateNodes(trees[[i]]))
        trees[[i]]$root.time <- oldRootAge + anchorTime 
      }
      return(trees)
    }
    #for a single tree
    if (class(trees) == "phylo"){
      oldRootAge <- max(dateNodes(trees))
      trees$root.time <- oldRootAge + anchorTime
    }
    return(trees)
}


get_Completeness <- function(p.file, p.stat = param){


# get Standard Deviation for completeness

# sort p.file for plotting
netDiv <- p.file[round(0.25*length(p.file[,1]) + 1):length(p.file[,1]),8:17]
turnover <- p.file[round(0.25*length(p.file[,1]) + 1):length(p.file[,1]),19:28]
relSamp <- p.file[round(0.25*length(p.file[,1]) + 1):length(p.file[,1]),30:39]

# convert FBD parameters to lambda, mu, and psi
meanMu <- calcMu(netDiv, turnover)
meanLambda <- calcLambda(netDiv, meanMu)
meanPsi <- calcPsi(meanMu, relSamp)

completeness <- meanPsi/ ((meanMu + meanPsi))

SD <- c()
for (i in 1:length(completeness[1,])){
	
	SD[i] <- sd(completeness[,i])
	
}

# read in the p.stats file
param <- read.table("Tip_dated_tree_analysis.nex.pstat", header = TRUE)

# sort 'param' into quantities to plot
netDiv <- param[5:14,]
upper.div <- param[5:14,5]
lower.div <- param[5:14,4]
turnover <- param[16:25,]
upper.t <- param[16:25,5]
lower.t <- param[16:25,4]
relSamp <- param[27:36,]
upper.s <- param[27:36,5]
lower.s <- param[27:36,4]

# convert FBD parameters to speciation, extinction, and sampling rates
meanMu <- calcMu(netDiv[,2], turnover[,2])
meanLambda <- calcLambda(netDiv[,2], meanMu)
meanPsi <- calcPsi(meanMu, relSamp[,2])

# if you want medians instead
meanMu <- calcMu(netDiv[,6], turnover[,6])
meanLambda <- calcLambda(netDiv[,6], meanMu)
meanPsi <- calcPsi(meanMu, relSamp[,6])

completeness <- meanPsi/ ((meanMu + meanPsi))

return(list("completeness" = completeness, "SD"= SD))
}


# plot completeness (Phi), where Phi = (Psi / (Psi + Mu))


Completeness_Through_Time <- function(midpoints, completeness, SD, up_to = 251){

plot(midpoints, completeness, xlim = c(485.4, up_to),ylim = c(-0.0005, 1), type = "n", axes = F, xlab = "", ylab = "")

rect(System[,1][-1], rep(0, 9), System[,1][-10], rep(100, 9), col = rep(c("grey95", "grey97")), border = NA)

segments(midpoints, completeness, midpoints, completeness+SD, col = "sienna1", lwd = 2)
segments(midpoints, completeness, midpoints, completeness-SD, col = "sienna1", lwd = 2)

lines(midpoints, completeness, xlim = c(485.4,0), lwd = 2, col = "sienna1")
points(midpoints, completeness, xlim = c(485.4,0), pch = 21, cex = 1.3, col = "black", bg = "sienna1")

axis(1, col = "grey75", at = c(485.4, 251), labels = FALSE)
axis(1, col = "grey75", at = seq(450 , up_to, by = -50))
axis(2, col = "grey75", line = 0.5, at = seq(0, 1, 0.2))

mtext("Age (Ma)", side = 1, line = 2)
mtext("Completeness", side = 2, line = 2.7)

}


# Get per-interval probability for sampling an ancestor-descendant pair using measure of completeness, equal to Phi^2

get_int_AD_pair_probs <- function(p.stat = param){
	
param <- read.table("Tip_dated_tree_analysis.nex.pstat", header = TRUE)

# sort 'param' into quantities to plot
netDiv <- param[5:14,]
upper.div <- param[5:14,5]
lower.div <- param[5:14,4]
turnover <- param[16:25,]
upper.t <- param[16:25,5]
lower.t <- param[16:25,4]
relSamp <- param[27:36,]
upper.s <- param[27:36,5]
lower.s <- param[27:36,4]

# convert FBD parameters to speciation, extinction, and sampling rates
meanMu <- calcMu(netDiv[,2], turnover[,2])
meanLambda <- calcLambda(netDiv[,2], meanMu)
meanPsi <- calcPsi(meanMu, relSamp[,2])

AD_pair <- (meanPsi/ ((meanMu + meanPsi)))^2; names(AD_pair) <- rownames(System); 

return(AD_pair)
}


# simulate a sample of brittle star morphologies from a single locality. Ideally, the probability distribuion is parameratized by the empirical data

sim.Sample <- function(N, m, sd){
	
# simulate size distributions for a fossil locality but make sure none of the trait values are less than zero (there's no such thing as negative size!)
		x <- 0
		while(x == 0){
			sim.param <- rnorm(n = N, mean = m, sd = sd)
			
			if(all(sim.param > 0)){
				x <- 1
			}
		}
	
		
	sim.mean <- mean(sim.param)
	sim.sd <- sd(sim.param)
	
	return(c(sim.mean, sim.sd))
}


# simulate a fossil record (i.e., time series) of brittle star morphologies. Here, Nsamples, M, and SD are vectors of # of samples, their means, and SD values, and need to be of the same length.

sim.fossil.Series <- function(Nsim, occ_times, Nsamples, M, SD){
	
	means <- matrix(data=NA, nrow = Nsim, ncol = length(occ_times))
	St.D <- matrix(data=NA, nrow = Nsim, ncol = length(occ_times))
	
	for (i in 1:Nsim){
	
		for (j in 1:length(occ_times)){
			N <- Nsamples[j]; m <- M[j]; sd <- SD[j]
			foo <- sim.Sample(N,m,sd)
			means[i,j] <- foo[1]
			St.D[i,j] <- foo[2]
		
		}
	
	}
	
	return(list("Mean" = means, "St.D" = St.D))
}

# testing stuff...
#Occs <- sort(runif(10, 1, 10))
#Nsamples <- round(runif(10, 1, 64))
#M <- rnorm(10, rnorm(1, 10, 1), 0.5)
#SD <- runif(10, 0.5, 1)

#k <- sim.fossil.Series(Nsim = 1, occ_times = Occs, Nsamples = Nsamples, M = M, SD = SD)
#plot(Occs, k$Mean, ylab = "Size (mm)", xlab = "Time (Ma)", ylim = c(min(k$Mean-k$St.D), max(k$Mean+k$St.D))); lines(Occs,k$Mean); points(Occs,k$Mean, cex = 1.2, pch = 19, col = "black"); segments(Occs, k$Mean, Occs, k$Mean+k$St.D);segments(Occs, k$Mean, Occs, k$Mean-k$St.D);


sim.Gotland.Series <- function(Nsim, occ_times, Nsamples, M, SD, bootstrap.M = FALSE, randomize.SD = FALSE){
	
	means <- matrix(data=NA, nrow = Nsim, ncol = length(occ_times))

	# simulate Nsim sequences
	for (i in 1:Nsim){
		
		#optional argument to bootstrap trait valuess over the time-series. Used for Monte Carlo tests of magnitude change and sustained directionality
		if(bootstrap.M == TRUE){
	
		M <- sample(M, length(M), replace = TRUE)
	
		}
		
		#optional argument to probablistically sample standard deviation
		if(randomize.SD == TRUE){
	
		SD <- runif(length(Nsamples), min(SD), max(SD))
	
		}
		
			# generate a value for each locality
			for (j in 1:length(occ_times)){
						
				N <- Nsamples[j]; m <- M[j]; sd <- SD[j]
				
				# probabilistically generate new trait values and store them
				foo <- sim.Sample(N,m,sd)
				means[i,j] <- foo[1]
			#	St.D[i,j] <- foo[2]
		
			}
	
		}
	
		
	return(list("sim.means" = means))
}


get_SE <- function(sim.means){
	
	St.E <- vector(mode = "numeric", length = length(occ_times))

	for (k in 1:length(occ_times)){
			
		# get standard error of the mean (i.e., St.D of the sampling distribution)
		St.E[k] <- sd(sim.means[,k])
			
	}
	
return(St.E)
	
}


magnitude.test <- function(sim.values){

IE <- c(); ME <- c(); LE <- c()

for (i in 1:Nsim){
	
	IE[i] <- sim.values$sim.means[i,1] - sim.values$sim.means[i,3]
	ME[i] <- sim.values$sim.means[i,4] - sim.values$sim.means[i,6]
	LE[i] <- sim.values$sim.means[i,7] - sim.values$sim.means[i,9]
	
}

Pr_IE <- length(which(IE <= M[1] - M[3])) / Nsim; 
Pr_ME <- length(which(IE <= M[4] - M[6])) / Nsim; 
Pr_LE <- length(which(IE <= M[7] - M[9])) / Nsim; 

print(paste("test statistic =", M[1] - M[3],", Probability of the observed magnitude of size decrease across the Ireviken Event = ", Pr_IE))
print(paste("test statistic =", M[4] - M[6],", Probability of the observed magnitude of size decrease across the Mulde Event = ", Pr_ME))
print(paste("test statistic =", M[7] - M[9],", Probability of the observed magnitude of size decrease across the Lau Event = ", Pr_LE))

}


# test for the directionality of sustained size decrease across Mulde and Lau Events

directional.test <- function(sim.values){

ME <- rep(0,Nsim); LE <- rep(0,Nsim)

	for (i in 1:Nsim){
	
		if (sim.values$sim.means[i,4] > sim.values$sim.means[i,5] && sim.values$sim.means[i,5] > sim.values$sim.means[i,6]){
		
		ME[i] <- 1
		
		}	
	
		if (sim.values$sim.means[i,7] > sim.values$sim.means[i,8] && sim.values$sim.means[i,8] > sim.values$sim.means[i,9]){
		
		LE[i] <- 1
		
		}		

	}

Pr_ME <- length(which(ME == 1)) / Nsim; 
Pr_LE <- length(which(LE == 1)) / Nsim; 

print(paste("Probability of sustained size decrease across the Mulde Event = ", Pr_ME))
print(paste("Probability of sustained size decrease across the Lau Event = ", Pr_LE))

}


