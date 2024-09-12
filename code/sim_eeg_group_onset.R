# Simulate ERP in Np participants, then estimate onsets using cluster sum.
# Function used in onsetsim_eeg_group_comp.Rmd
sim_eeg_group_onset <- function(ronset, Nt, Np, ...){
  
  out <- vector(mode = "numeric", length = Np)
  
  cond1 <- matrix(0, nrow = Nt, ncol = Nf)
  cond2 <- matrix(0, nrow = Nt, ncol = Nf)
  
  for(P in 1:Np){ # participants
    
    ponset <- sample(ronset, 1) # get random onset
    st <- which(Xf==ponset) # find starting point
    temp2 <- c(rep(0, st-2), erp, rep(0, Nf-st-length(erp)+2)) # pad vector
    
    for(T in 1:Nt){
      cond2[T,] <- temp2 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
      cond1[T,] <- temp1 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
    }
    
    # t-tests
    ori.t2 <- vector(mode = "numeric", length = Nf)
    for(F in 1:Nf){
      ori.t2[F] <- t.test(cond1[,F], cond2[,F])$statistic^2
    }
    
    # Make permutation table of t values 
    perm.t2 <- permtdist(cond1, cond2, Nt, Nf, nboot = nboot)^2
    perm.th <- apply(perm.t2, 2, quantile, probs = 1-aath)
    
    perm.pvals <- vector(mode = "numeric", length = Nf)
    for(F in 1:Nf){
      perm.pvals[F] <- (sum(perm.t2[,F] >= ori.t2[F]) + 1) / (nboot + 1)
    }
    
    # cluster-sum statistics -----
    cmap <- cluster.make(perm.pvals <= aath)
    perm.max.sums <- vector(mode = "numeric", length = nboot)
    for(B in 1:nboot){
      # threshold permutation t2 values and form clusters
      perm.cmap <- cluster.make(perm.t2[B,] <= perm.th)  
      perm.max.sums[B] <- max(cluster.sum(values = perm.t2[B,], cmap = perm.cmap))
    }
    # cluster sum threshold
    cs.th <- quantile(perm.max.sums, probs = 1-aath)
    # cluster test
    cs.test <- cluster.test(values = ori.t2, cmap = cmap, cs.th)
    out[P] <- find_onset(cs.test, Xf)
  }
  out # return vector of onsets
}

# Simulate ERP in Np participants, then estimate onsets.
# Use change point detection instead of cluster sum.
sim_eeg_group_onset_cp <- function(ronset, Nt, Np, ...){
  
  out <- vector(mode = "numeric", length = Np)
  
  cond1 <- matrix(0, nrow = Nt, ncol = Nf)
  cond2 <- matrix(0, nrow = Nt, ncol = Nf)
  
  for(P in 1:Np){ # participants
    
    ponset <- sample(ronset, 1) # get random onset
    st <- which(Xf==ponset) # find starting point
    temp2 <- c(rep(0, st-2), erp, rep(0, Nf-st-length(erp)+2)) # pad vector
    
    for(T in 1:Nt){
      cond2[T,] <- temp2 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
      cond1[T,] <- temp1 + eeg_noise(frames = Nf, srate = srate, outvar = outvar, meanpower)
    }
    
    # t-tests
    ori.t2 <- vector(mode = "numeric", length = Nf)
    for(F in 1:Nf){
      ori.t2[F] <- t.test(cond1[,F], cond2[,F])$statistic^2
    }
    
    # fit change point model
    res <- cpt.meanvar(ori.t2, method = "BinSeg", Q=2, param.estimates = FALSE)
    out[P] <- Xf[res@cpts[1]]
    
  }
  out # return vector of onsets
}
