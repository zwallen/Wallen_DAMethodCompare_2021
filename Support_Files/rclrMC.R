# Robust Centered Log-Ratio Transformation with Matrix Completion
# using OptSpace Algorithm (rclrMC)
#
# R implementation of the robust centered log-ratio transformation
# with matrix completion using OptSpace algorithm detailed in the manuscript
# Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A, 
# Knight R, Zengler K. 2019. A novel sparse compositional 
# technique reveals microbial perturbations. mSystems 4:e00016-19. 
#
# Parameters:
# phy_obj -- A phyloseq object that contains taxonomic
#            abundances (raw counts, i.e. not transformed or normalized
#            for inter-sample differences in sequencing depth) that can be
#            accessed using phyloseq's otu_table() function.
# ropt -- Set NA to guess the rank, or a positive integer as a pre-defined rank.
#         See OptSpace() function from ROptSpace package for more details.
# niter -- Maximum number of iterations allowed.
#          See OptSpace() function from ROptSpace package for more details.
# tol -- Stopping criterion for reconstruction in Frobenius norm.
#        See OptSpace() function from ROptSpace package for more details.
# seed -- A number to pass to set.seed() function, so transformations can
#         be reproduced.

rclrMC = function(phy_obj, var, control, ropt = NA, niter = 5, tol = 1e-5, seed=123){
  
  # Load required libraries
  if (!require(phyloseq)){install.packages('phyloseq')}
  suppressMessages(library(phyloseq))
  if (!require(ROptSpace)){install.packages('ROptSpace')}
  suppressMessages(library(ROptSpace))
  
  # Check that input is a phyloseq object that has taxonomic abundances
  if (tryCatch(dim(otu_table(phy_obj))[1], error=function(x){return("error")}) == "error"){
    stop("ERROR: input is not a phyloseq object with taxa abundances accessible using otu_table()")
  }
  
  # Check that there are no negative values
  if (any(otu_table(phy_obj) < 0)){
    stop("ERROR: input contains negative values")
  }
  
  # Check that no undefined values exist (Inf or -Inf)
  if (any(otu_table(phy_obj) == Inf) || any(otu_table(phy_obj) == -Inf)){
    stop("ERROR: input contains Inf and/or -Inf values")
  }
  
  # Check for missing data
  if (any(is.na(otu_table(phy_obj)))){
    stop("ERROR: input contains NAs or NaNs")
  }
  
  # Check to make sure zeros are present in data
  if (!any(otu_table(phy_obj) == 0)){
    stop("ERROR: input does not contain any zeros, therefore, standard clr transformation should be applied")
  }
  
  # Make sure samples are rows and taxa are columns
  if (taxa_are_rows(phy_obj)){
    phy_obj <- phyloseq(t(otu_table(phy_obj)))
  }
  
  # Check for samples that contain only zeros
  if (any(rowSums(otu_table(phy_obj)) == 0)){
    stop("ERROR: input contains samples with all zeros")
  }

  # Convert counts to sample centered data (via total sum scaling)
  phy_obj <- transform_sample_counts(phy_obj, function(x){x / sum(x)})
  
  # Ensure all centered data sums to 1
  if (any(round(rowSums(otu_table(phy_obj)),0) != 1)){
    stop("ERROR: total sum scaling of counts did not result in samples adding to 1")
  }
  
  # Take log of sample centered data
  phy_obj <- transform_sample_counts(phy_obj, function(x){log(x)})
  
  # Mask non-finite values
  otu_table(phy_obj)[!is.finite(otu_table(phy_obj))] <- NA
  
  # Calculate centered log-ratios centering data on geometric mean
  phy_obj <- transform_sample_counts(phy_obj, function(x){x - mean(na.omit(x))})
  
  # Perform OptSpace algorithm for matrix completion
  set.seed(seed)
  opt_out <- OptSpace(matrix(otu_table(phy_obj), nrow = dim(otu_table(phy_obj))[1], ncol = dim(otu_table(phy_obj))[2]),
                      ropt = ropt, niter = niter, tol = tol, showprogress = FALSE)
  
  # Replace missing data with imputed data
  idf <- data.frame(opt_out$X %*% opt_out$S %*% t(opt_out$Y))
  otu_table(phy_obj)[is.na(otu_table(phy_obj))] <- idf[is.na(otu_table(phy_obj))]
  
  # Return phyloseq object with transformed values
  return(phy_obj)
}
  
