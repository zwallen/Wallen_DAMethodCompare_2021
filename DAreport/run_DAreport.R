#!/usr/bin/env Rscript
#
# This program will run DAreport on the command line using the same parameters
# as the function itself. The usual pdf reports will be produced, but the results
# for each method and other results will be outputted into a multi-tab excel file
# titled DAreport_results.xlsx. If listing more than one method to skip make sure
# it is a comma separated list with no spaces.
#
# Usage:
#    ./run_DAreport ps.rds group_variable_name prev_cut method1,method2,etc 1234

###################################################################################
# DAreport (Differential abundance report)                                        #
#                                                                                 #
# R function to run various DA methods compared in Wallen 2021 BMC Bioinformatics #
# on a feature abundance table and generate a report showing the distribution     #
# of concordances between DA method results (DAreport_concordance.pdf) and what   #
# differentially abundant features are being detected in relation to their mean   #
# relative abundance and fold change between two groups of interest               #
# (DAreport_MRA_vs_FC.pdf). Results of DA methods and other data used in the      #
# will be returned in a list.                                                     #
#                                                                                 #
# Note: certain modifications have been made since the original BMC Bioinformatics#
# comparison study. Modifications include:                                        #
#           - Log TSS is now being performed by normalizing using TSS first, then #
#             log transforming the normalized values after replacing zeros with   #
#             a pseudocount of half the minimum TSS normalized value as this way  #
#             is more intuitive and is the way it is done in popular methods such #
#             as MaAsLin2.                                                        #
#           - Robust CLR with matrix completion is now using default parameters   #
#             for the OptSpace() function.                                        #
#           - ANCOM v2 R code has been updated to ANCOM v2.1 R code.              #
#             See FrederickHuangLin/ANCOM/scripts/ancom_v2.1.R.                   #
#           - The strategy and functions used for GLM NBZI has changed. Now using #
#             glm.nb() function to perform negative binomial GLM and zeroinfl()   #
#             function to perform zero-inflated negative binomial GLM. This helps #
#             with computation time. NA pvalues are given a 1.                    #
#           - Wilcoxon rank-sum test is now being used instead of Kruskal-Wallis  #
#             as this would be the more appropriate test to use with only two     #
#             groups of interest.                                                 #
#                                                                                 #
# Note: make sure LEfSe scripts are accessible via your PATH variable. I found it #
# easiest to download LEfSe via making its own conda environment and running      #
# conda install -c biobakery lefse to download it and its dependencies.           #
#                                                                                 #
# Parameters:                                                                     #
#            ps - Input phyloseq object that contains the observed, untransformed #
#                 or normalized (to total sequence depth) abundances of features  #
#                 (taxa, gene families, pathways, etc.). Should have an otu_table #
#                 component with the abundances, and a sample_data component with #
#                 a variable defining the two groups of interest. This variable   #
#                 should be coded as 1 (for treatment/experimental group) and 0   #
#                 (for control group).                                            #
#     group.var - The variable name in the sample_data component of the phyloseq  #
#                 object that contains data on which samples belong to which of   #
#                 the two groups of interest.                                     #
#      prev.cut - The prevalence cutoff for features to be included in the DA     #
#                 testing. Should be a value between 0 and 1, i.e. the proportion #
#                 of samples a feature should be present in to include. Default is#
#                 1 to include all features.                                      #
#          skip - Vector of name(s) of DA method(s) that you would like to skip.  #
#          seed - Numeric value to pass to the set.seed() function in order to    #
#                 to make reports consistent.                                     #
###################################################################################

DAreport <- function(ps, group.var, prev.cut=1, skip=NULL, seed=1234){

  # load required libraries
  suppressWarnings(suppressMessages(library(phyloseq)))
  suppressWarnings(suppressMessages(library(dplyr)))
  suppressWarnings(suppressMessages(library(ROptSpace)))
  suppressWarnings(suppressMessages(library(ALDEx2)))
  suppressWarnings(suppressMessages(library(ANCOMBC)))
  suppressWarnings(suppressMessages(library(compositions)))
  suppressWarnings(suppressMessages(library(nlme)))
  suppressWarnings(suppressMessages(library(DESeq2)))
  suppressWarnings(suppressMessages(library(pscl)))
  suppressWarnings(suppressMessages(library(MASS)))
  suppressWarnings(suppressMessages(library(samr)))
  suppressWarnings(suppressMessages(library(edgeR)))
  suppressWarnings(suppressMessages(library(baySeq)))
  suppressWarnings(suppressMessages(library(metagenomeSeq)))
  suppressWarnings(suppressMessages(library(limma)))
  suppressWarnings(suppressMessages(library(ggplot2)))
  suppressWarnings(suppressMessages(library(reshape2)))
  suppressWarnings(suppressMessages(library(ggh4x)))
  suppressWarnings(suppressMessages(library(ggpubr)))
  suppressWarnings(suppressMessages(library(grid)))
  suppressWarnings(suppressMessages(library(scales)))

  ###############################################################
  ########### define ANCOM v2.1 functions #######################
  # copied from FrederickHuangLin/ANCOM/scripts/ancom_v2.1.R
  # and added here so source file does not have to be pointed to
  ###############################################################
  # Data Pre-Processing
  feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL,
                                       out_cut = 0.05, zero_cut = 0.90, lib_cut, neg_lb){
    feature_table = data.frame(feature_table, check.names = FALSE)
    meta_data = data.frame(meta_data, check.names = FALSE)
    # Drop unused levels
    meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
    # Match sample IDs between metadata and feature table
    sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
    feature_table = feature_table[, sample_ID]
    meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]

    # 1. Identify outliers within each taxon
    if (!is.null(group_var)) {
      group = meta_data[, group_var]
      z = feature_table + 1 # Add pseudo-count (1)
      f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
      f_fit = lm(f ~ group)
      e = rep(0, length(f)); e[!is.na(group)] = residuals(f_fit)
      y = t(t(z) - e)

      outlier_check = function(x){
        # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
        mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
        sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
        pi = 0.75
        n = length(x)
        epsilon = 100
        tol = 1e-5
        score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
        while (epsilon > tol) {
          grp1_ind = (score >= 1)
          mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
          sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
          sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
          pi_new = sum(grp1_ind)/n

          para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
          if(any(is.na(para))) break

          score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
            ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))

          epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 +
                           (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
          mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
        }

        if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
          if(pi < out_cut){
            out_ind = grp1_ind
          }else if(pi > 1 - out_cut){
            out_ind = (!grp1_ind)
          }else{
            out_ind = rep(FALSE, n)
          }
        }else{
          out_ind = rep(FALSE, n)
        }
        return(out_ind)
      }
      out_ind = matrix(FALSE, nrow = nrow(feature_table), ncol = ncol(feature_table))
      out_ind[, !is.na(group)] = t(apply(y, 1, function(i)
        unlist(tapply(i, group, function(j) outlier_check(j)))))

      feature_table[out_ind] = NA
    }

    # 2. Discard taxa with zeros  >=  zero_cut
    zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
    taxa_del = which(zero_prop >= zero_cut)
    if(length(taxa_del) > 0){
      feature_table = feature_table[- taxa_del, ]
    }

    # 3. Discard samples with library size < lib_cut
    lib_size = colSums(feature_table, na.rm = T)
    if(any(lib_size < lib_cut)){
      subj_del = which(lib_size < lib_cut)
      feature_table = feature_table[, - subj_del]
      meta_data = meta_data[- subj_del, ]
    }

    # 4. Identify taxa with structure zeros
    if (!is.null(group_var)) {
      group = factor(meta_data[, group_var])
      present_table = as.matrix(feature_table)
      present_table[is.na(present_table)] = 0
      present_table[present_table != 0] = 1

      p_hat = t(apply(present_table, 1, function(x)
        unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
      samp_size = t(apply(feature_table, 1, function(x)
        unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
      p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)

      struc_zero = (p_hat == 0) * 1
      # Whether we need to classify a taxon into structural zero by its negative lower bound?
      if(neg_lb) struc_zero[p_hat_lo <= 0] = 1

      # Entries considered to be structural zeros are set to be 0s
      struc_ind = struc_zero[, group]
      feature_table = feature_table * (1 - struc_ind)

      colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
    }else{
      struc_zero = NULL
    }

    # 5. Return results
    res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero)
    return(res)
  }

  # ANCOM main function
  ANCOM = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH",
                   alpha = 0.05, adj_formula = NULL, rand_formula = NULL, ...){
    # OTU table transformation:
    # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
    if (!is.null(struc_zero)) {
      num_struc_zero = apply(struc_zero, 1, sum)
      comp_table = feature_table[num_struc_zero == 0, ]
    }else{
      comp_table = feature_table
    }
    comp_table = log(as.matrix(comp_table) + 1)
    n_taxa = dim(comp_table)[1]
    taxa_id = rownames(comp_table)
    n_samp = dim(comp_table)[2]

    # Determine the type of statistical test and its formula.
    if (is.null(rand_formula) & is.null(adj_formula)) {
      # Basic model
      # Whether the main variable of interest has two levels or more?
      if (length(unique(meta_data%>%pull(main_var))) == 2) {
        # Two levels: Wilcoxon rank-sum test
        tfun = stats::wilcox.test
      } else{
        # More than two levels: Wilcoxon test
        tfun = stats::kruskal.test
      }
      # Formula
      tformula = formula(paste("x ~", main_var, sep = " "))
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      # Model: ANOVA
      tfun = stats::aov
      # Formula
      tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
    }else if (!is.null(rand_formula)) {
      # Model: Mixed-effects model
      tfun = nlme::lme
      # Formula
      if (is.null(adj_formula)) {
        # Random intercept model
        tformula = formula(paste("x ~", main_var))
      }else {
        # Random coefficients/slope model
        tformula = formula(paste("x ~", main_var, "+", adj_formula))
      }
    }

    # Calculate the p-value for each pairwise comparison of taxa.
    p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
    colnames(p_data) = taxa_id
    rownames(p_data) = taxa_id
    for (i in 1:(n_taxa - 1)) {
      # Loop through each taxon.
      # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
      # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted
      # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
      # the first (i - 1) columns were calculated by previous iterations, and
      # the i^th column contains all zeros.
      alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ])
      # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
      # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
      alr_data = alr_data[, - (1:i), drop = FALSE]
      n_lr = dim(alr_data)[2] # number of log-ratios (lr)
      alr_data = cbind(alr_data, meta_data) # merge with the metadata

      # P-values
      if (is.null(rand_formula) & is.null(adj_formula)) {
        p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
          suppressWarnings(tfun(tformula,
                                data = data.frame(x, alr_data,
                                                  check.names = FALSE))$p.value)
        }
        )
      }else if (is.null(rand_formula) & !is.null(adj_formula)) {
        p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
          fit = tfun(tformula,
                     data = data.frame(x, alr_data, check.names = FALSE),
                     na.action = na.omit)
          summary(fit)[[1]][main_var, "Pr(>F)"]
        }
        )
      }else if (!is.null(rand_formula)) {
        p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
          fit = tfun(fixed = tformula,
                     data = data.frame(x, alr_data, check.names = FALSE),
                     random = formula(rand_formula),
                     na.action = na.omit, ...)
          anova(fit)[main_var, "p-value"]
        }
        )
      }
    }
    # Complete the p-value matrix.
    # What we got from above iterations is a lower triangle matrix of p-values.
    p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
    diag(p_data) = 1 # let p-values on diagonal equal to 1
    p_data[is.na(p_data)] = 1 # let p-values of NA equal to 1

    # Multiple comparisons correction.
    q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))

    # Calculate the W statistic of ANCOM.
    # For each taxon, count the number of q-values < alpha.
    W = apply(q_data, 2, function(x) sum(x < alpha))

    # Organize outputs
    out_comp = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
    # Declare a taxon to be differentially abundant based on the quantile of W statistic.
    # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
    out_comp = out_comp%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
                                 detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
                                 detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
                                 detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))

    # Feature with structural zeros are automatically declared to be differentially abundant
    if (!is.null(struc_zero)){
      out = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE,
                       detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE,
                       row.names = NULL, check.names = FALSE)
      out[match(taxa_id, out$taxa_id), ] = out_comp
    }else{
      out = out_comp
    }

    # Draw volcano plot
    # Calculate clr
    clr_table = apply(feature_table, 2, clr)
    # Calculate clr mean difference
    eff_size = apply(clr_table, 1, function(y)
      lm(y ~ x, data = data.frame(y = y,
                                  x = meta_data %>% pull(main_var),
                                  check.names = FALSE))$coef[-1])

    if (is.matrix(eff_size)){
      # Data frame for the figure
      dat_fig = data.frame(taxa_id = out$taxa_id, t(eff_size), y = out$W, check.names = FALSE) %>%
        mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No"))) %>%
        gather(key = group, value = x, rownames(eff_size))
      # Replcace "x" to the name of covariate
      dat_fig$group = sapply(dat_fig$group, function(x) gsub("x", paste0(main_var, " = "), x))
      # Replace Inf by (n_taxa - 1) for structural zeros
      dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

      fig = ggplot(data = dat_fig) + aes(x = x, y = y) +
        geom_point(aes(color = zero_ind)) +
        facet_wrap(~ group) +
        labs(x = "CLR mean difference", y = "W statistic") +
        scale_color_discrete(name = "Structural zero", drop = FALSE) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "top",
              strip.background = element_rect(fill = "white"))
      fig
    } else{
      # Data frame for the figure
      dat_fig = data.frame(taxa_id = out$taxa_id, x = eff_size, y = out$W) %>%
        mutate(zero_ind = factor(ifelse(is.infinite(y), "Yes", "No"), levels = c("Yes", "No")))
      # Replace Inf by (n_taxa - 1) for structural zeros
      dat_fig$y = replace(dat_fig$y, is.infinite(dat_fig$y), n_taxa - 1)

      fig = ggplot(data = dat_fig) + aes(x = x, y = y) +
        geom_point(aes(color = zero_ind)) +
        labs(x = "CLR mean difference", y = "W statistic") +
        scale_color_discrete(name = "Structural zero", drop = FALSE) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.position = "top")
      fig
    }

    res = list(out = out, fig = fig)
    return(res)
  }
  ########### end ANCOM functions ###########
  ###########################################

  ##################################################
  ########### define robust CLR function ###########
  # Robust Centered Log-Ratio Transformation with Matrix Completion
  # using OptSpace algorithm (rclrMC)
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
  #    ropt -- Set NA to guess the rank, or a positive integer as a pre-defined rank.
  #            See OptSpace() function from ROptSpace package for more details.
  #   niter -- Maximum number of iterations allowed.
  #            See OptSpace() function from ROptSpace package for more details.
  #     tol -- Stopping criterion for reconstruction in Frobenius norm.
  #            See OptSpace() function from ROptSpace package for more details.
  #    seed -- A number to pass to set.seed() function, so transformations can
  #            be reproduced.

  rclrMC = function(phy_obj, ropt = NA, niter = 50, tol = 1e-6, seed=1234){

    # Check that input is a phyloseq object that has feature abundances
    if (tryCatch(dim(otu_table(phy_obj))[1], error=function(x){return("error")}) == "error"){
      stop("input is not a phyloseq object with feature abundances accessible using otu_table()")
    }

    # Check that there are no negative values
    if (any(otu_table(phy_obj) < 0)){
      stop("input contains negative values")
    }

    # Check that no undefined values exist (Inf or -Inf)
    if (any(otu_table(phy_obj) == Inf) || any(otu_table(phy_obj) == -Inf)){
      stop("input contains Inf and/or -Inf values")
    }

    # Check for missing data
    if (any(is.na(otu_table(phy_obj)))){
      stop("input contains NAs or NaNs")
    }

    # Make sure samples are rows and features are columns
    if (taxa_are_rows(phy_obj)){
      phy_obj <- phyloseq(t(otu_table(phy_obj)))
    }

    # Check for samples that contain only zeros
    if (any(rowSums(otu_table(phy_obj)) == 0)){
      stop("input contains samples with all zeros")
    }

    # CLR transform using geometric mean of non-zero components
    phy_obj <- transform_sample_counts(phy_obj, function(x){log(x)-mean(log(x[x>0]))})

    # Make non-finite values NA
    otu_table(phy_obj)[!is.finite(otu_table(phy_obj))] <- NA

    # Perform OptSpace algorithm for matrix completion
    set.seed(seed)
    opt_out <- OptSpace(matrix(otu_table(phy_obj), nrow = nrow(otu_table(phy_obj)), ncol = ncol(otu_table(phy_obj))),
                        ropt = ropt, niter = niter, tol = tol, showprogress = FALSE)

    # Replace missing data with imputed data
    idf <- data.frame(opt_out$X %*% opt_out$S %*% t(opt_out$Y))
    otu_table(phy_obj)[is.na(otu_table(phy_obj))] <- idf[is.na(otu_table(phy_obj))]

    # Return phyloseq object with transformed values
    return(phy_obj)
  }
  ########### end robust CLR function code  ###########
  #####################################################

  #################################################
  ########### perform some data checks  ###########

  # check to make input is valid phyloseq object with data
  if (tryCatch(nrow(otu_table(ps)), error=function(x){return('error')}) == 'error'){
    stop('phyloseq object supplied is either not a phyloseq object, or has no otu_table data')
  }

  # check that there are no negative values
  if (any(otu_table(ps) < 0)){
    stop("input abundances contains negative values")
  }

  # check that no undefined values exist (Inf or -Inf)
  if (any(otu_table(ps) == Inf) || any(otu_table(ps) == -Inf)){
    stop("input abundances contains Inf and/or -Inf values")
  }

  # check for missing data
  if (any(is.na(otu_table(ps)))){
    stop("input abundances contains NAs or NaNs")
  }

  # check for samples that contain only zeros
  if (any(rowSums(otu_table(ps)) == 0)){
    stop("input abundances contains samples with all zeros")
  }

  # check to make sure group variable is present in phyloseq object
  suppressWarnings(
  if (tryCatch(sample_variables(ps), error=function(x){return('error')}) == 'error'){
    stop('no sample variables detected in supplied phyloseq object')
  }else if (!(group.var %in% sample_variables(ps))){
    stop('given group variable name was not found in column names of sample data')
  })

  # check to make sure group variable has only two levels and they are 0 and 1
  group.var.levels <- levels(factor(as.character(data.frame(sample_data(ps)[,group.var])[,1])))
  if (length(group.var.levels) != 2){
    stop('group variable does not have 2 levels, only grouping variables with 2 levels is supported at this time')
  }
  if (!("0" %in% group.var.levels)){
    stop('0 should be one of the two levels in the specified group level')
  }
  if (!("1" %in% group.var.levels)){
    stop('1 should be one of the two levels in the specified group level')
  }
  if (length(grep("0|1", group.var.levels, invert = T)) != 0){
    stop('group variable should only contain values of 0 and 1')
  }

  # make sure samples are rows and features are columns
  if (taxa_are_rows(ps)){
    ps <- phyloseq(t(otu_table(ps)))
  }

  # make sure feature names are R acceptable
  taxa_names(ps) <- make.names(taxa_names(ps))

  # make group variable a factor variable
  sample_data(ps)[,group.var] <- factor(as.character(data.frame(sample_data(ps)[,group.var])[,1]))

  ########### end data checks  ###########
  ########################################

  ######################################################################
  ##### perform needed transformations and calculate effect sizes  #####

  cat('\n', '### Performing needed data transformations ###', '\n')

  # tss
  if (length(grep('_tss', skip))==0){
    cat('\n', 'total sum scaling (TSS)...', '\n')
    ps.tss <- transform_sample_counts(ps, function(x){x/sum(x)})
  }

  # log tss
  if (length(grep('_tss', skip))==0){
    cat('\n', 'log TSS...', '\n')
    log.trans <- function(x) {
      y <- replace(x, x == 0, min(x[x>0]) / 2)
      return(log(y))
    }
    ps.log.tss <- ps.tss
    otu_table(ps.log.tss) <- otu_table(log.trans(data.frame(otu_table(ps.tss))), taxa_are_rows=FALSE)
  }

  # clr
  if (length(grep('_clr', skip))==0){
    cat('\n', 'centered log-ratio (CLR)...', '\n')
    ps.clr <- transform_sample_counts(ps, function(x){log(x+1)-mean(log(x+1))})
  }

  # rclr
  if (length(grep('_rclr', skip))==0){
    cat('\n', 'robust CLR with matrix completion...', '\n')
    ps.rclr <- rclrMC(ps, seed = seed)
  }

  ########### end tranformations  ###########
  ###########################################

  ################################################
  ##### perform requested feature filtering  #####

  cat('\n', '### Performing requested data filtering ###', '\n')

  ps <- filter_taxa(ps, function(x){sum(x > 0) > (1-prev.cut)*length(x)}, TRUE)
  if (length(grep('_tss', skip))==0){
    ps.tss <- prune_taxa(taxa_names(ps.tss)[taxa_names(ps.tss) %in% taxa_names(ps)], ps.tss)
    ps.log.tss <- prune_taxa(taxa_names(ps.log.tss)[taxa_names(ps.log.tss) %in% taxa_names(ps)], ps.log.tss)
  }
  if (length(grep('_clr', skip))==0){
    ps.clr <- prune_taxa(taxa_names(ps.clr)[taxa_names(ps.clr) %in% taxa_names(ps)], ps.clr)
  }
  if (length(grep('_rclr', skip))==0){
    ps.rclr <- prune_taxa(taxa_names(ps.rclr)[taxa_names(ps.rclr) %in% taxa_names(ps)], ps.rclr)
  }

  cat('\n', 'Features found in less than', (1-prev.cut)*nsamples(ps), 'samples have been removed', '\n')
  cat(' Features remaining for differential abundance testing:',ntaxa(ps), '\n')

  ########### end filtering  ###########
  ######################################

  ############################################################
  ###### begin running of differential abundance tests  ######

  cat('\n', '### Running differential abundance methods ###', '\n')

  # initialize data structures to capture outputs to be used later
  result.list <- list()
  conc.input <- data.frame(row.names = taxa_names(ps))

  #### ALDEx2 ####
  if (!('aldex2' %in% skip)){
    cat('\n', 'Running ALDEx2...', '\n')
    # set the comparison groups
    conds <- as.vector(data.frame(sample_data(ps)[,group.var])[,1])
    # perfrom ALDEx2 analysis
    set.seed(seed)
    aldex2 <- aldex(t(data.frame(round(otu_table(ps), 0))),
                    conds,
                    mc.samples = 1000,
                    test = "t",
                    effect =FALSE,
                    include.sample.summary = FALSE,
                    denom = "all",
                    verbose = TRUE)
    # extract results
    aldex2 <- data.frame(Feature=rownames(aldex2), aldex2, row.names = NULL, check.names = FALSE)
    # make sure all features accounted for
    if (nrow(aldex2) != length(taxa_names(ps))){
      miss.feat <- taxa_names(ps)[!(make.names(taxa_names(ps)) %in% aldex2$Feature)]
      aldex2 <- rbind(aldex2, data.frame(Feature=miss.feat, we.ep=rep(1, length(miss.feat)),
                                         we.eBH=rep(1, length(miss.feat)), wi.ep=rep(1, length(miss.feat)),
                                         wi.eBH=rep(1, length(miss.feat))))
    }
    # append
    aldex2 <- aldex2[match(rownames(conc.input), aldex2$Feature),]
    result.list <- c(result.list, aldex2=list(aldex2))
    conc.input$`ALDEx2 t-test`[aldex2$we.eBH < 0.05] <- 1
    conc.input$`ALDEx2 Wilcoxon`[aldex2$wi.eBH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### ANCOM ####
  if (!('ancom' %in% skip)){
    cat('\n', 'Running ANCOM...', '\n')
    # extract feature count data from phyloseq object
    feat.df <- data.frame(t(otu_table(ps)), check.names = FALSE)
    # extract sample metadata from phyloseq object
    meta.df <- data.frame(Sample.ID=sample_names(ps), sample_data(ps)[,group.var], check.names = FALSE)
    meta.df[,group.var] <- factor(as.character(meta.df[,group.var]))
    # run data pre-processing function
    pre.process <- feature_table_pre_process(feat.df, meta.df, sample_var = "Sample.ID",
                                             group_var = NULL, out_cut = 0, zero_cut = 1,
                                             lib_cut = 0, neg_lb = FALSE)
    # perform ANCOM
    ancom.out <- ANCOM(feature_table = pre.process$feature_table,
                       meta_data = pre.process$meta_data,
                       struc_zero = pre.process$structure_zeros,
                       main_var = group.var,
                       p_adj_method = "BH",
                       alpha = 0.05,
                       adj_formula = NULL,
                       rand_formula = NULL)
    # extract results
    ancom <- ancom.out$out
    colnames(ancom)[1] <- "Feature"
    # append
    ancom <- ancom[match(rownames(conc.input), ancom$Feature),]
    result.list <- c(result.list, ancom=list(ancom))
    conc.input$ANCOM[ancom$detected_0.8 == TRUE] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### ANCOM-BC ####
  if (!('ancombc' %in% skip)){
    cat('\n', 'Running ANCOM-BC...', '\n')
    # perform ANCOM-BC
    ancom.out <- ancombc(phyloseq = ps,
                         formula = group.var,
                         p_adj_method = "BH",
                         zero_cut = 1,
                         lib_cut = 0,
                         group = NULL,
                         struc_zero = FALSE,
                         neg_lb = FALSE,
                         tol = 1e-05,
                         max_iter = 100,
                         conserve = FALSE,
                         alpha = 0.05,
                         global = FALSE)
    # extract results
    ancombc <- data.frame(Feature=rownames(ancom.out$res$beta), Beta=ancom.out$res$beta[,1],
                          SE=ancom.out$res$se[,1], P=ancom.out$res$p_val[,1], FDR_BH=ancom.out$res$q_val[,1],
                          row.names = NULL, check.names = FALSE)
    # append
    ancombc <- ancombc[match(rownames(conc.input), ancombc$Feature),]
    result.list <- c(result.list, ancombc=list(ancombc))
    conc.input$`ANCOM-BC`[ancombc$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### DESeq2 nbinomWaldTest ####
  if (!('deseq2_wald' %in% skip)){
    cat('\n', 'Running DESeq2 nbinomWaldTest...', '\n')
    # convert phyloseq object to DESeqDataSet object specifying design formula with only case/control status as predictor
    dds <- phyloseq_to_deseq2(ps, formula(paste("~",group.var)))
    # calculate size factors using modified geometric means that can handle zeros and estimate normalization factors
    dds <- estimateSizeFactors(dds, type="poscounts")
    # perform differential abundance analysis
    dds <- DESeq(dds)
    # extract results
    deseq2_wald <- results(dds, cooksCutoff = FALSE)
    deseq2_wald <- data.frame(Feature=rownames(deseq2_wald), deseq2_wald, row.names = NULL, check.names = FALSE)
    # append
    deseq2_wald <- deseq2_wald[match(rownames(conc.input), deseq2_wald$Feature),]
    result.list <- c(result.list, deseq2_wald=list(deseq2_wald))
    conc.input$`DESeq2 nbinomWaldTest`[deseq2_wald$padj < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### negative binomial GLM ####
  if (!('glm_nbzi' %in% skip)){
    cat('\n', 'Running negative binomial GLM...', '\n')
    # calculate total sequence count to use as offset variable in model
    total_seq_count <- as.vector(rowSums(otu_table(ps)))
    # loop negative binomial GLM model for all features taking results of zero-inflated
    # model if AIC is smaller than model without zero-inflation, else taking results
    # of model without zero-inflation
    glm_nbzi <- data.frame()
    for (i in 1:length(taxa_names(ps))){
      if (sum(round(as.vector(otu_table(ps)[,i]),0)) > 0){
        if (min(otu_table(ps)[,i]) > 0){
          mod_zero <- NULL
        }else{
          mod_zero <- suppressWarnings(suppressMessages(
                       zeroinfl(round(as.vector(otu_table(ps)[,i]),0) ~ data.frame(sample_data(ps)[,group.var])[,1] +
                                                                        offset(log(total_seq_count)),
                                dist = "negbin", link = "logit")))
        }

        mod <- tryCatch(glm.nb(round(as.vector(otu_table(ps)[,i]), 0) ~ data.frame(sample_data(ps)[,group.var])[,1] +
                                                                        offset(log(total_seq_count))),
                        error=function(e){"err"})

        if (length(mod) == 1){
          beta <- summary(mod_zero)$coefficients$count[2,1]
          se <- summary(mod_zero)$coefficients$count[2,2]
          pval <- summary(mod_zero)$coefficients$count[2,4]
          zero <- TRUE
        }else if (is.null(mod_zero)){
          beta <- summary(mod)$coefficients[2,1]
          se <- summary(mod)$coefficients[2,2]
          pval <- summary(mod)$coefficients[2,4]
          zero <- FALSE
        }else if (is.na(summary(mod_zero)$coefficients$count[2,4])){
          beta <- summary(mod)$coefficients[2,1]
          se <- summary(mod)$coefficients[2,2]
          pval <- summary(mod)$coefficients[2,4]
          zero <- FALSE
        }else if (AIC(mod_zero) < AIC(mod)){
          beta <- summary(mod_zero)$coefficients$count[2,1]
          se <- summary(mod_zero)$coefficients$count[2,2]
          pval <- summary(mod_zero)$coefficients$count[2,4]
          zero <- TRUE
        }else{
          beta <- summary(mod)$coefficients[2,1]
          se <- summary(mod)$coefficients[2,2]
          pval <- summary(mod)$coefficients[2,4]
          zero <- FALSE
        }
      }else{
        beta <- NA
        se <- NA
        pval <- 1
        zero <- NA
      }
      glm_nbzi <- rbind(glm_nbzi, data.frame(Feature = taxa_names(ps)[i],
                                             Zero_Infl = zero,
                                             Beta = beta,
                                             SE = se,
                                             P = pval))
    }
    # make any NA p-values 1
    glm_nbzi$P[is.na(glm_nbzi$P)] <- 1
    # perform FDR correction for pvalues
    glm_nbzi$FDR_BH <- p.adjust(glm_nbzi$P, method = 'BH')
    # append
    glm_nbzi <- glm_nbzi[match(rownames(conc.input), glm_nbzi$Feature),]
    result.list <- c(result.list, glm_nbzi=list(glm_nbzi))
    conc.input$`GLM NBZI`[glm_nbzi$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### LEfSe ####
  if (!('lefse' %in% skip)){
    cat('\n', 'Running LEfSe...', '\n')
    # extract feature count data
    lefse.df <- data.frame(otu_table(ps.tss), row.names=NULL)
    # add column for sample and group designations
    lefse.df <- cbind(data.frame(sample=sample_names(ps.tss),
                                 group=data.frame(sample_data(ps.tss)[,group.var])[,1]), lefse.df)
    # write out LEfSe input data
    write.table(lefse.df, "LEfSe_input.txt", quote = F, sep="\t", row.names = F)
    ### run LEfSe workflow ###
    # requires for LEfSe scripts and their dependencies to be in your PATH,
    # which I found is easiest by making a conda environment
    # for LEfSe with Python 2 and downloading through conda
    # run python script for formatting
    system("lefse_format_input.py LEfSe_input.txt LEfSe_input.in -f c -c 2 -s -1 -u 1 -o 1000000")
    # run LEfSe using default parameters
    system("lefse_run.py LEfSe_input.in LEfSe.txt")
    # run LEfSe using parameters that will output all LDA scores and pvalues
    system("lefse_run.py LEfSe_input.in LEfSe_full.txt -a 1.0 -w 1.0 -l 0.0")
    # clean up input files
    system("rm LEfSe_input.txt")
    system("rm LEfSe_input.in")
    ##########################
    # read in results
    lefse <- read.table("LEfSe.txt", sep="\t")
    lefse.full <- read.table("LEfSe_full.txt")
    colnames(lefse) <- c("Feature", "LogOfHighestClassAvg", "HighestClass", "LDAScore", "P")
    # add FDR q-values in LEfSe results and further filter for those that are FDR < 0.05
    lefse$FDR_BH <- p.adjust(lefse.full[,5], method="BH")
    lefse[lefse$FDR_BH > 0.05,c(5,6)] <- "-"
    lefse[lefse$P == "-",6] <- "-"
    lefse[is.na(lefse$LDAScore),6] <- "-"
    lefse[lefse[,6] == "-",3] <- "-"
    lefse[lefse[,6] == "-",4] <- "-"
    # clean up remaining files
    system("rm LEfSe.txt")
    system("rm LEfSe_full.txt")
    # append
    lefse <- lefse[match(rownames(conc.input), lefse$Feature),]
    result.list <- c(result.list, lefse=list(lefse))
    conc.input$LEfSe[lefse$FDR_BH != "-"] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### SAMseq ####
  if (!('samseq' %in% skip)){
    cat('\n', 'Running SAMseq...', '\n')
    # recode group variable with 1 and 2 (required by SAMseq)
    groups <- recode(data.frame(sample_data(ps)[,group.var])[,1], "1"=2L, "0"=1L)
    # normalize data using SAM method and round to whole integer (required by SAMseq)
    norm.data <- round(samr.norm.data(t(otu_table(ps))), 0)
    # run SAMseq analysis
    set.seed(seed)
    sam.results <- SAMseq(norm.data, groups, resp.type = "Two class unpaired", geneid = taxa_names(ps), fdr.output = 1)
    # create results table
    samseq <- data.frame(Feature = rownames(sam.results$samr.obj$x), FoldChange = sam.results$samr.obj$foldchange)
    sig.gene.table <- rbind(sam.results$siggenes.table$genes.up, sam.results$siggenes.table$genes.lo)[,-1]
    colnames(sig.gene.table)[1] <- "Feature"
    samseq <- merge(data.frame(samseq), data.frame(sig.gene.table, check.names = F), by = "Feature", all = TRUE)
    # append
    samseq <- samseq[match(rownames(conc.input), samseq$Feature),]
    result.list <- c(result.list, samseq=list(samseq))
    conc.input$SAMseq[as.numeric(samseq$`q-value(%)`) < 5 & !is.na(samseq$`q-value(%)`)] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### baySeq ####
  if (!('bayseq' %in% skip)){
    cat('\n', 'Running baySeq...', '\n')
    # calculate total sequence count to use as size factor for baySeq
    total_seq_count <- as.vector(rowSums(otu_table(ps)))
    # combine data into countData object and supply library size
    replicates <- data.frame(sample_data(ps)[,group.var])[,1]
    groups <- list(NDE = as.character(rep(0, dim(otu_table(ps))[1])),
                   DE = data.frame(sample_data(ps)[,group.var])[,1])
    count.data <- round(t(data.frame(otu_table(ps))),0)
    cd <- new("countData", data=count.data, groups=groups, replicates=replicates)
    cd@annotation  <- data.frame(name = taxa_names(ps))
    libsizes(cd) <- total_seq_count
    # estimate prior parameters using default negative binomial distribution
    cd <- getPriors.NB(cd, cl = NULL)
    # find posterior likelihood of each count or paired count belonging to one of the defined models
    cd <- getLikelihoods(cd, cl = NULL)
    # extract results
    bayseq <- topCounts(cd, group = "DE", number = Inf)
    bayseq <- data.frame(Feature=bayseq$name, bayseq[,c("likes","DE","FDR.DE","FWER.DE")],
                         row.names = NULL, check.names = FALSE)
    # append
    bayseq <- bayseq[match(rownames(conc.input), bayseq$Feature),]
    result.list <- c(result.list, bayseq=list(bayseq))
    conc.input$baySeq[bayseq$FDR.DE < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### edgeR exactTest with RLE ####
  if (!('edger_exactTest_rle' %in% skip)){
    cat('\n', 'Running edgeR exactTest with RLE...', '\n')
    # convert phyloseq object to DGEList object specifying groups to be tested
    dge <- DGEList(counts=t(as(otu_table(ps), "matrix"))+1,
                   group=data.frame(sample_data(ps)[,group.var])[,1],
                   genes=taxa_names(ps))
    # calculate normalization factors and tagwise dispersions
    dge <- edgeR::calcNormFactors(dge, method="RLE")
    dge <- estimateTagwiseDisp(estimateCommonDisp(dge))
    # perform differential abundance analysis
    et.results <- exactTest(dge)
    # FDR correct Pvalues and extract results
    edger_exactTest_rle <- data.frame(topTags(et.results, n = nrow(dge$table), adjust.method = "BH", sort.by = "none")$table,
                                      row.names = NULL, check.names = FALSE)
    colnames(edger_exactTest_rle)[1] <- "Feature"
    # append
    edger_exactTest_rle <- edger_exactTest_rle[match(rownames(conc.input), edger_exactTest_rle$Feature),]
    result.list <- c(result.list, edger_exactTest_rle=list(edger_exactTest_rle))
    conc.input$`edgeR exactTest RLE`[edger_exactTest_rle$FDR < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### edgeR exactTest with TMM ####
  if (!('edger_exactTest_tmm' %in% skip)){
    cat('\n', 'Running edgeR exactTest with TMM...', '\n')
    # convert phyloseq object to DGEList object specifying groups to be tested
    dge <- DGEList(counts=t(as(otu_table(ps), "matrix")),
                   group=data.frame(sample_data(ps)[,group.var])[,1],
                   genes=taxa_names(ps))
    # calculate normalization factors and tagwise dispersions
    dge <- edgeR::calcNormFactors(dge)
    dge <- estimateTagwiseDisp(estimateCommonDisp(dge))
    # perform differential abundance analysis
    et.results <- exactTest(dge)
    # FDR correct Pvalues and extract results
    edger_exactTest_tmm <- data.frame(topTags(et.results, n = nrow(dge$table), adjust.method = "BH", sort.by = "none")$table,
                            row.names = NULL, check.names = FALSE)
    colnames(edger_exactTest_tmm)[1] <- "Feature"
    # append
    edger_exactTest_tmm <- edger_exactTest_tmm[match(rownames(conc.input), edger_exactTest_tmm$Feature),]
    result.list <- c(result.list, edger_exactTest_tmm=list(edger_exactTest_tmm))
    conc.input$`edgeR exactTest TMM`[edger_exactTest_tmm$FDR < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### fitFeatureModel ####
  if (!('fitfeat' %in% skip)){
    cat('\n', 'Running fitFeatureModel...', '\n')
    # convert phyloseq object to MRexperiment object
    mre <- phyloseq_to_metagenomeSeq(ps)
    # perform Cummulative Sum Scaling transformation on data
    mre <- cumNorm(mre)
    # perform differential abundance analysis
    mod <- with(pData(mre), model.matrix(formula(paste("~",group.var))))
    fit.results <- fitFeatureModel(mre, mod = mod)
    # FDR correct Pvalues and extract results
    fitfeat <- MRfulltable(fit.results, number = nrow(assayData(mre)$counts), adjustMethod = "BH")
    fitfeat <- data.frame(Feature=rownames(fitfeat), fitfeat, row.names = NULL, check.names = FALSE)
    # append
    fitfeat <- fitfeat[match(rownames(conc.input), fitfeat$Feature),]
    result.list <- c(result.list, fitfeat=list(fitfeat))
    conc.input$fitFeatureModel[fitfeat$adjPvalues < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### fitZIG ####
  if (!('fitzig' %in% skip)){
    cat('\n', 'Running fitZIG...', '\n')
    # convert phyloseq object to MRexperiment object
    mre <- phyloseq_to_metagenomeSeq(ps)
    # perform Cummulative Sum Scaling transformation on data
    mre <- cumNorm(mre)
    # perform differential abundance analysis
    mod <- with(pData(mre), model.matrix(formula(paste("~",group.var))))
    fit.results <- fitZig(mre, mod = mod, useCSSoffset = T, control = zigControl(dfMethod="default"))
    # FDR correct Pvalues and extract results
    fitzig <- MRfulltable(fit.results, number = nrow(assayData(mre)$counts), adjustMethod = "BH")
    fitzig <- data.frame(Feature=rownames(fitzig), fitzig, row.names = NULL, check.names = FALSE)
    # append
    fitzig <- fitzig[match(rownames(conc.input), fitzig$Feature),]
    result.list <- c(result.list, fitzig=list(fitzig))
    conc.input$fitZIG[fitzig$adjPvalues < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### limma voom ####
  if (!('limma_voom' %in% skip)){
    cat('\n', 'Running limma voom...', '\n')
    # convert phyloseq object to DGEList object specifying groups to be tested
    dge <- DGEList(counts=t(as(otu_table(ps), "matrix")),
                   group=data.frame(sample_data(ps)[,group.var])[,1],
                   genes=taxa_names(ps))
    # create model matrix
    mod <- with(data.frame(sample_data(ps), check.names = FALSE), model.matrix(formula(paste("~",group.var))))
    # use voom to transform to log2-CPM and calculate associated weights to ready data for linear modeling with limma
    dge.voom <- voom(dge, design=mod, plot = FALSE)
    # perform linear modeling and calculate empirical Bayes statistics for differential expression
    fit <- lmFit(dge.voom, design=mod)
    fit <- eBayes(fit)
    # FDR correct Pvalues and extract results
    limma_voom <- data.frame(topTable(fit, coef=paste(group.var,"1",sep=""), n = Inf, sort.by = "none", adjust.method="BH"),
                                 row.names = NULL, check.names = FALSE)
    colnames(limma_voom)[1] <- "Feature"
    # append
    limma_voom <- limma_voom[match(rownames(conc.input), limma_voom$Feature),]
    result.list <- c(result.list, limma_voom=list(limma_voom))
    conc.input$`limma-voom`[limma_voom$adj.P.Val < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### t-test ####
  if (!('t_tss' %in% skip)){
    cat('\n', 'Running t-test with TSS...', '\n')
    # Loop t-test for all features
    t_tss <- data.frame()
    for (i in 1:length(taxa_names(ps.log.tss))){
      res <- t.test(as.vector(otu_table(ps.log.tss)[,i]) ~ data.frame(sample_data(ps.log.tss)[,group.var])[,1])
      t_tss <- rbind(t_tss, data.frame(Feature = taxa_names(ps.log.tss)[i],
                                       Group_1_mean = res$estimate[1],
                                       Group_0_mean = res$estimate[2],
                                       Mean_Ratio = res$estimate[1]/res$estimate[2],
                                       P = res$p.value, row.names = NULL))
    }
    # perform FDR correction for pvalues
    t_tss$FDR_BH <- p.adjust(t_tss$P, method = 'BH')
    # append
    t_tss <- t_tss[match(rownames(conc.input), t_tss$Feature),]
    result.list <- c(result.list, t_tss=list(t_tss))
    conc.input$`t-test TSS`[t_tss$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }
  if (!('t_clr' %in% skip)){
    cat('\n', 'Running t-test with CLR...', '\n')
    # Loop t-test for all features
    t_clr <- data.frame()
    for (i in 1:length(taxa_names(ps.clr))){
      res <- t.test(as.vector(otu_table(ps.clr)[,i]) ~ data.frame(sample_data(ps.clr)[,group.var])[,1])
      t_clr <- rbind(t_clr, data.frame(Feature = taxa_names(ps.clr)[i],
                                       Group_1_mean = res$estimate[1],
                                       Group_0_mean = res$estimate[2],
                                       Mean_Ratio = res$estimate[1]/res$estimate[2],
                                       P = res$p.value, row.names = NULL))
    }
    # perform FDR correction for pvalues
    t_clr$FDR_BH <- p.adjust(t_clr$P, method = 'BH')
    # append
    t_clr <- t_clr[match(rownames(conc.input), t_clr$Feature),]
    result.list <- c(result.list, t_clr=list(t_clr))
    conc.input$`t-test CLR`[t_clr$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }
  if (!('t_rclr' %in% skip)){
    cat('\n', 'Running t-test with robust CLR...', '\n')
    # Loop t-test for all features
    t_rclr <- data.frame()
    for (i in 1:length(taxa_names(ps.rclr))){
      res <- t.test(as.vector(otu_table(ps.rclr)[,i]) ~ data.frame(sample_data(ps.rclr)[,group.var])[,1])
      t_rclr <- rbind(t_rclr, data.frame(Feature = taxa_names(ps.rclr)[i],
                                         Group_1_mean = res$estimate[1],
                                         Group_0_mean = res$estimate[2],
                                         Mean_Ratio = res$estimate[1]/res$estimate[2],
                                         P = res$p.value, row.names = NULL))
    }
    # perform FDR correction for pvalues
    t_rclr$FDR_BH <- p.adjust(t_rclr$P, method = 'BH')
    # append
    t_rclr <- t_rclr[match(rownames(conc.input), t_rclr$Feature),]
    result.list <- c(result.list, t_rclr=list(t_rclr))
    conc.input$`t-test rCLR`[t_rclr$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### Wicoxon rank sum test ####
  if (!('w_tss' %in% skip)){
    cat('\n', 'Running Wilcoxon with TSS...', '\n')
    # Loop Wilcoxon for all features
    w_tss <- data.frame()
    for (i in 1:length(taxa_names(ps.tss))){
      group_1_mean <- mean(otu_table(ps.tss)[sample_data(ps.tss)[,group.var]=="1",][,i])
      group_0_mean <- mean(otu_table(ps.tss)[sample_data(ps.tss)[,group.var]=="0",][,i])
      res <- wilcox.test(as.vector(otu_table(ps.tss)[,i]) ~ data.frame(sample_data(ps.tss)[,group.var])[,1])
      w_tss <- rbind(w_tss, data.frame(Feature = taxa_names(ps.tss)[i],
                                         Group_1_mean = group_1_mean,
                                         Group_0_mean = group_0_mean,
                                         Mean_Ratio = group_1_mean/group_0_mean,
                                         P = res$p.value, row.names = NULL))
    }
    # perform FDR correction for pvalues
    w_tss$FDR_BH <- p.adjust(w_tss$P, method = 'BH')
    # append
    w_tss <- w_tss[match(rownames(conc.input), w_tss$Feature),]
    result.list <- c(result.list, w_tss=list(w_tss))
    conc.input$`Wilcoxon TSS`[w_tss$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }
  if (!('w_clr' %in% skip)){
    cat('\n', 'Running Wilcoxon with CLR...', '\n')
    # Loop Wilcoxon for all features
    w_clr <- data.frame()
    for (i in 1:length(taxa_names(ps.clr))){
      group_1_mean <- mean(otu_table(ps.clr)[sample_data(ps.clr)[,group.var]=="1",][,i])
      group_0_mean <- mean(otu_table(ps.clr)[sample_data(ps.clr)[,group.var]=="0",][,i])
      res <- wilcox.test(as.vector(otu_table(ps.clr)[,i]) ~ data.frame(sample_data(ps.clr)[,group.var])[,1])
      w_clr <- rbind(w_clr, data.frame(Feature = taxa_names(ps.clr)[i],
                                         Group_1_mean = group_1_mean,
                                         Group_0_mean = group_0_mean,
                                         Mean_Ratio = group_1_mean/group_0_mean,
                                         P = res$p.value, row.names = NULL))
    }
    # perform FDR correction for pvalues
    w_clr$FDR_BH <- p.adjust(w_clr$P, method = 'BH')
    # append
    w_clr <- w_clr[match(rownames(conc.input), w_clr$Feature),]
    result.list <- c(result.list, w_clr=list(w_clr))
    conc.input$`Wilcoxon CLR`[w_clr$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }
  if (!('w_rclr' %in% skip)){
    cat('\n', 'Running Wilcoxon with robust CLR...', '\n')
    # Loop Wilcoxon for all features
    w_rclr <- data.frame()
    for (i in 1:length(taxa_names(ps.rclr))){
      group_1_mean <- mean(otu_table(ps.rclr)[sample_data(ps.rclr)[,group.var]=="1",][,i])
      group_0_mean <- mean(otu_table(ps.rclr)[sample_data(ps.rclr)[,group.var]=="0",][,i])
      res <- wilcox.test(as.vector(otu_table(ps.rclr)[,i]) ~ data.frame(sample_data(ps.rclr)[,group.var])[,1])
      w_rclr <- rbind(w_rclr, data.frame(Feature = taxa_names(ps.rclr)[i],
                                           Group_1_mean = group_1_mean,
                                           Group_0_mean = group_0_mean,
                                           Mean_Ratio = group_1_mean/group_0_mean,
                                           P = res$p.value, row.names = NULL))
    }
    # perform FDR correction for pvalues
    w_rclr$FDR_BH <- p.adjust(w_rclr$P, method = 'BH')
    # append
    w_rclr <- w_rclr[match(rownames(conc.input), w_rclr$Feature),]
    result.list <- c(result.list, w_rclr=list(w_rclr))
    conc.input$`Wilcoxon rCLR`[w_rclr$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  #### GLM ####
  if (!('glm_tss' %in% skip)){
    cat('\n', 'Running GLM with TSS...', '\n')
    # Loop GLM for all features
    glm_tss <- data.frame()
    for (i in 1:length(taxa_names(ps.log.tss))){
      group_1_mean <- mean(otu_table(ps.log.tss)[sample_data(ps.log.tss)[,group.var]=="1",][,i])
      group_0_mean <- mean(otu_table(ps.log.tss)[sample_data(ps.log.tss)[,group.var]=="0",][,i])
      res <- glm(as.vector(otu_table(ps.log.tss)[,i]) ~ data.frame(sample_data(ps.log.tss)[,group.var])[,1])
      glm_tss <- rbind(glm_tss, data.frame(Feature = taxa_names(ps.log.tss)[i],
                                           Group_1_mean = group_1_mean,
                                           Group_0_mean = group_0_mean,
                                           Mean_Ratio = group_1_mean/group_0_mean,
                                           Beta = summary(res)$coefficients[2,1],
                                           SE = summary(res)$coefficients[2,2],
                                           P = summary(res)$coefficients[2,4],
                                           row.names = NULL))
    }
    # perform FDR correction for pvalues
    glm_tss$FDR_BH <- p.adjust(glm_tss$P, method = 'BH')
    # append
    glm_tss <- glm_tss[match(rownames(conc.input), glm_tss$Feature),]
    result.list <- c(result.list, glm_tss=list(glm_tss))
    conc.input$`GLM TSS`[glm_tss$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }
  if (!('glm_clr' %in% skip)){
    cat('\n', 'Running GLM with CLR...', '\n')
    # Loop GLM for all features
    glm_clr <- data.frame()
    for (i in 1:length(taxa_names(ps.clr))){
      group_1_mean <- mean(otu_table(ps.clr)[sample_data(ps.clr)[,group.var]=="1",][,i])
      group_0_mean <- mean(otu_table(ps.clr)[sample_data(ps.clr)[,group.var]=="0",][,i])
      res <- glm(as.vector(otu_table(ps.clr)[,i]) ~ data.frame(sample_data(ps.clr)[,group.var])[,1])
      glm_clr <- rbind(glm_clr, data.frame(Feature = taxa_names(ps.clr)[i],
                                           Group_1_mean = group_1_mean,
                                           Group_0_mean = group_0_mean,
                                           Mean_Ratio = group_1_mean/group_0_mean,
                                           Beta = summary(res)$coefficients[2,1],
                                           SE = summary(res)$coefficients[2,2],
                                           P = summary(res)$coefficients[2,4],
                                           row.names = NULL))
    }
    # perform FDR correction for pvalues
    glm_clr$FDR_BH <- p.adjust(glm_clr$P, method = 'BH')
    # append
    glm_clr <- glm_clr[match(rownames(conc.input), glm_clr$Feature),]
    result.list <- c(result.list, glm_clr=list(glm_clr))
    conc.input$`GLM CLR`[glm_clr$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }
  if (!('glm_rclr' %in% skip)){
    cat('\n', 'Running GLM with robust CLR...', '\n')
    # Loop GLM for all features
    glm_rclr <- data.frame()
    for (i in 1:length(taxa_names(ps.rclr))){
      group_1_mean <- mean(otu_table(ps.rclr)[sample_data(ps.rclr)[,group.var]=="1",][,i])
      group_0_mean <- mean(otu_table(ps.rclr)[sample_data(ps.rclr)[,group.var]=="0",][,i])
      res <- glm(as.vector(otu_table(ps.rclr)[,i]) ~ data.frame(sample_data(ps.rclr)[,group.var])[,1])
      glm_rclr <- rbind(glm_rclr, data.frame(Feature = taxa_names(ps.rclr)[i],
                                             Group_1_mean = group_1_mean,
                                             Group_0_mean = group_0_mean,
                                             Mean_Ratio = group_1_mean/group_0_mean,
                                             Beta = summary(res)$coefficients[2,1],
                                             SE = summary(res)$coefficients[2,2],
                                             P = summary(res)$coefficients[2,4],
                                             row.names = NULL))
    }
    # perform FDR correction for pvalues
    glm_rclr$FDR_BH <- p.adjust(glm_rclr$P, method = 'BH')
    # append
    glm_rclr <- glm_rclr[match(rownames(conc.input), glm_rclr$Feature),]
    result.list <- c(result.list, glm_rclr=list(glm_rclr))
    conc.input$`GLM rCLR`[glm_rclr$FDR_BH < 0.05] <- 1
    conc.input[is.na(conc.input)] <- 0
  }

  ########### end differential abundance tests  ###########
  #########################################################

  ################################################
  ########### begin report generation  ###########

  cat('\n', '### Generating reports ###', '\n')

  #### concordance portion of report ####
  cat('\n', 'Generating concordance based report...', '\n')

  # calculate pairwise concordances for method results
  conc.df <- apply(conc.input, 2,
                   function(x){apply(conc.input, 2,
                                     function(y){sum(x + y == 2)/sum(x + y > 0)})
                   }
  )
  diag(conc.df) <- NA

  # order methods by mean concordance
  means <- sort(apply(conc.df, 2, function(x){mean(na.omit(x))}))
  conc.df <- conc.df[match(names(means), rownames(conc.df)), match(names(means), colnames(conc.df))]

  # calculate relative abundances and fold changes
  group_1_mra <- apply(otu_table(ps.tss)[sample_data(ps.tss)[,group.var]=="1",], 2, mean)
  group_0_mra <- apply(otu_table(ps.tss)[sample_data(ps.tss)[,group.var]=="0",], 2, mean)
  conc.input$MRA <- (group_1_mra+group_0_mra)/2
  conc.input$FC <- group_1_mra/group_0_mra

  # calculate proportion of differential abundance signatures
  # detected by each method out of the total features tested
  prop.det <- round(colSums(conc.input[,1:(ncol(conc.input)-2)])/nrow(conc.input), 2)

  # melt data
  suppressWarnings(suppressMessages(melt.conc <- melt(conc.df)))

  # add calculated proportions
  melt.conc$prop <- NA
  for (i in 1:length(prop.det)){
    melt.conc$prop[melt.conc[,2] == names(prop.det)[i]] <- prop.det[i]
  }

  # double data making one half data for boxplots, and the other half for bar plot
  melt.conc <- rbind(data.frame(melt.conc, plot="1"), data.frame(melt.conc, plot="2"))
  melt.conc$value[melt.conc$plot == 2] <- NA

  # add mean concordance to plot for boxplots of concordances
  melt.conc$mean <- NA
  melt.conc$mean[melt.conc$plot == 1] <- mean(na.omit(melt.conc$value))

  # add color for plotting
  col.pal <- rainbow(n = length(levels(melt.conc[,2])))
  set.seed(seed)
  names(col.pal) <- sample(levels(melt.conc[,2]))

  melt.conc$color <- NA
  for (i in 1:nrow(melt.conc)){
    melt.conc[i,"color"] <- col.pal[as.character(melt.conc[i,2])]
  }

  # get correlation labels for dot plot
  conc.prop.cor <- cor.test(melt.conc$prop, melt.conc$value)
  cor.lab <- data.frame(label = paste("Pearson's r [95% CI] = ",
                                      round(conc.prop.cor$estimate,2), " [",
                                      round(conc.prop.cor$conf.int[1],2), "; ",
                                      round(conc.prop.cor$conf.int[2],2), "]",", ",
                                      "P = ",signif(conc.prop.cor$p.value,0),sep=""))

  # generate plots
  g1.1 <- ggplot(data = melt.conc[melt.conc$plot==1,], aes(x = Var2, y = value)) +
    geom_boxplot(fill="grey", notch=F, outlier.size=-1) +
    geom_dotplot(inherit.aes=F, data = melt.conc[melt.conc$plot==1,], aes(x = Var2, y = value), binaxis = "y",
                 stackdir = "center", binwidth = 0.03, position = position_dodge(0.75), dotsize = 0.3) +
    stat_summary(fun=mean, geom="point", shape=1, size=1, color="red") +
    geom_hline(data=melt.conc[melt.conc$plot==1,], aes(yintercept=mean), color="red", linetype="solid", size=0.25) +
    geom_bar(inherit.aes=F, data=unique(melt.conc[melt.conc$plot==2,2:ncol(melt.conc)]), aes(x = Var2, y = prop, fill=Var2), stat="identity") +
    facet_nested(plot ~ ., labeller=labeller(plot=c("1"="Concordances", "2"="Proportion DA"))) +
    scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
    scale_x_discrete(position="top") +
    scale_fill_manual(values=unique(melt.conc[melt.conc$plot==1,]$color)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 25, hjust = 0, size = 5), axis.text.x.top = element_text(vjust = 0.05)) +
    theme(axis.text.y = element_text(size = 6)) +
    theme(strip.text = element_text(margin = ggplot2::margin(0,0.01,0,0.01, "in"), size = 7)) +
    theme(plot.margin = unit(c(-0.75,1.5,1.25,-0.5), "lines")) +
    guides(fill="none") +
    labs(x="", y="")

  suppressWarnings(
  g1.2 <- ggplot(data = melt.conc[melt.conc$plot==1,], aes(x = prop, y = value)) +
    stat_smooth(method="lm", n=nrow(melt.conc[melt.conc$plot==1,]), color="black", alpha=0.8) +
    geom_point(aes(color = Var2), size=1) +
    scale_color_manual(values=unique(melt.conc[melt.conc$plot==1,]$color)) +
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    geom_text(data=cor.lab, aes(x=0.5, y=1.02, label=label), size=2) +
    scale_x_continuous(position="top") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 6), axis.text.x = element_text(size = 6)) +
    theme(axis.title.y = element_text(size = 6), axis.text.y = element_text(size = 6)) +
    theme(strip.text = element_text(margin = ggplot2::margin(0,0.01,0,0.01, "in"), size = 7)) +
    theme(plot.margin = unit(c(0.7,0.5,1.25,-0.5), "lines")) +
    theme(legend.position="right", legend.title=element_blank(), legend.text = element_text(size=6),
          legend.spacing.x = unit(0.01, "in"),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10)) +
    guides(color=guide_legend(reverse=T, ncol=1, keyheight = 0.5, keywidth = 0.5, override.aes = list(size=2))) +
    labs(x="Proportion DA", y="Concordances")
  )

  suppressMessages(suppressWarnings(g1 <- ggarrange(g1.1, g1.2, ncol=2, labels=NULL)))

  # add calculated concordances and proportions of DA features to method results list
  result.list <- c(result.list, concordances=list(conc.df), prop_DA=list(data.frame(Method=names(prop.det), Proportion_DA=prop.det)))

  #### differentially abundant features as a function of mean relative abundance ####
  #### and effect size (fold change)                                         ####
  cat('\n', 'Generating report on DA features as a function of MRA and fold change...', '\n')

  # prep binary feature by method matrix data for plotting
  suppressWarnings(suppressMessages(mra.plot <- melt(conc.input[,1:(ncol(conc.input)-2)])))
  mra.plot$MRA <- rep(conc.input$MRA, nrow(conc.df))
  mra.plot$FC <- rep(conc.input$FC, nrow(conc.df))
  mra.plot <- mra.plot[is.finite(mra.plot$FC),]
  mra.plot <- mra.plot[mra.plot$FC > 0,]

  # Order based on previous figure
  mra.plot$variable <- factor(mra.plot$variable, levels=rev(levels(melt.conc$Var2)))

  # add median mean relative abundance for plotting
  mra.plot$median <- NA
  mra.plot$median <- median(na.omit(mra.plot$MRA))

  # add variable that tags what MRA is above or below median
  mra.plot$prev <- ifelse(mra.plot$MRA > mra.plot$median, 1, 0)

  # add variable that tags what fold change is within or outside a log2 fold change of 0.4 (~1.3x change)
  mra.plot$fc_area <- ifelse(log2(mra.plot$FC) < -0.4 | log2(mra.plot$FC) > 0.4, 1, 0)

  # perform Fisher exact test for enrichment of certain feature as DA signatures based on MRA and log2(FC)
  mra.plot.lab <- data.frame()
  for (i in 1:length(levels(mra.plot$variable))){
    curr.method <- levels(mra.plot$variable)[i]
    sub.data <- mra.plot[mra.plot$variable == curr.method,]
    f.out.prev <- tryCatch(fisher.test(table(sub.data$prev, sub.data$value)), error=function(e){'err'})
    f.out.fc <- tryCatch(fisher.test(table(sub.data$fc_area, sub.data$value)), error=function(e){'err'})
    if (f.out.prev == "err"){
      f.out.prev <- list(estimate=NA, p.value=NA)
    }
    if (f.out.fc == "err"){
      f.out.fc <- list(estimate=NA, p.value=NA)
    }
    mra.plot.lab <- rbind(mra.plot.lab, data.frame(variable=curr.method,
                                                   prev_lab=paste("MRA: OR=",round(f.out.prev$estimate, 2),", P=",signif(f.out.prev$p.value, 0),sep=""),
                                                   fc_lab=paste("FC: OR=",round(f.out.fc$estimate, 2),", P=",signif(f.out.fc$p.value, 0),sep="")))
  }
  mra.plot.lab$variable <- factor(mra.plot.lab$variable, levels=levels(mra.plot$variable))
  result.list <- c(result.list, mra_fc_enrichment_tests=list(data.frame(Method=mra.plot.lab$variable,
                                                                        MRA_fisher_results=mra.plot.lab$prev_lab,
                                                                        FC_fisher_results=mra.plot.lab$fc_lab)))

  # generate plots
  g2 <- ggplot(data = mra.plot[order(mra.plot$value),], aes(x = log2(FC), y = MRA, color=as.character(value))) +
    geom_point(size=0.25) +
    geom_vline(data=mra.plot[order(mra.plot$value),], aes(xintercept=0.4), color="red", linetype="dashed", size=0.25) +
    geom_vline(data=mra.plot[order(mra.plot$value),], aes(xintercept=-0.4), color="red", linetype="dashed", size=0.25) +
    geom_hline(data=mra.plot[order(mra.plot$value),], aes(yintercept=median), color="red", linetype="dashed", size=0.25) +
    geom_text(inherit.aes=F, data=mra.plot.lab, aes(x=max(log2(mra.plot$FC)), y=1, label=prev_lab), size=1.5, position = position_dodge(width = 5)) +
    geom_text(inherit.aes=F, data=mra.plot.lab, aes(x=max(log2(mra.plot$FC)), y=0.1, label=fc_lab), size=1.5, position = position_dodge(width = 5)) +
    facet_grid(variable ~ .) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values=c("#999999", "#000000")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 5), axis.text.x = element_text(size = 4)) +
    theme(axis.title.y = element_text(size = 5), axis.text.y = element_text(size = 4)) +
    theme(strip.text = element_text(margin = ggplot2::margin(0,0.01,0,0.01, "in"), size = 3.5)) +
    theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
    theme(panel.spacing = unit(0.1, "lines")) +
    guides(color="none") +
    labs(x="log2 fold change in group 1", y="Mean relative abundances (log scale)")

  ########### end generating reports  ###########
  ###############################################

  ######################################
  ########### write reports  ###########

  pdf('DAreport_concordance.pdf', height=6, width=12)
  print(g1)
  dev.off()

  pdf('DAreport_MRA_vs_FC.pdf', height=nrow(conc.df), width=3)
  print(g2)
  dev.off()

  cat('\n', '### Done ###', '\n')

  return(result.list)
}

# load needed packages
suppressMessages(suppressWarnings(library(openxlsx)))

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
ps <- readRDS(args[1])
group.var <- args[2]
prev.cut <- args[3]
skip <- strsplit(args[4], ",")[[1]]
seed <- args[5]

# run DAreport
da.report <- suppressWarnings(DAreport(ps, group.var, prev.cut, skip, seed))

# print result list to excel
wb <- createWorkbook()
for (i in 1:length(da.report)){
addWorksheet(wb, names(da.report)[i])
writeData(wb, names(da.report)[i], da.report[[i]], keepNA=T)
}
saveWorkbook(wb, "DAreport_results.xlsx", overwrite=T)
