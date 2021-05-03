### Perform differential abundance testing using negative binomial GLM with or without zero-inflation

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(glmmTMB); cat("Running glmmTMB package version:", "\n"); packageVersion("glmmTMB")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Calculate total sequence number to use as offset variable in model
sample_data(ps)$total_seq_count <- rowSums(otu_table(ps))

# Dummy code case/control status
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)

# Loop negative binomial GLM model for all taxa taking results of zero-inflated
# model if AIC is smaller than model without zero-inflation, else taking results
# of model without zero-inflation
results <- data.frame()
for (i in 1:length(taxa_names(ps))){

  asv <- taxa_names(ps)[i]
  
  mod <- glmmTMB(as.vector(otu_table(ps)[,i]) ~ sample_data(ps)$case_control + offset(log(sample_data(ps)$total_seq_count)),
                   family = "nbinom2", ziformula = ~0)
  mod_zero <- glmmTMB(as.vector(otu_table(ps)[,i]) ~ sample_data(ps)$case_control + offset(log(sample_data(ps)$total_seq_count)),
                        family = "nbinom2", ziformula = ~1)

  if (is.na(AIC(mod_zero))){
    aic <- AIC(mod)
    beta <- summary(mod)$coefficients$cond[2,1]
    se <- summary(mod)$coefficients$cond[2,2]
    pval <- summary(mod)$coefficients$cond[2,4]
    zero <- FALSE
  }else if (AIC(mod_zero) < AIC(mod)){
    aic <- AIC(mod_zero)
    beta <- summary(mod_zero)$coefficients$cond[2,1]
    se <- summary(mod_zero)$coefficients$cond[2,2]
    pval <- summary(mod_zero)$coefficients$cond[2,4]
    zero <- TRUE
  }else{
    aic <- AIC(mod)
    beta <- summary(mod)$coefficients$cond[2,1]
    se <- summary(mod)$coefficients$cond[2,2]
    pval <- summary(mod)$coefficients$cond[2,4]
    zero <- FALSE
  }
  
  results <- rbind(results, data.frame(Representative_ASV = asv, 
                                       Zero_Infl = zero,
                                       Beta = beta,
				       SE = se,
				       P = pval))
}

# Perform FDR correction for pvalues
results$FDR_BH <- p.adjust(results$P, method = 'BH')

# Add taxa designations to table
results <- inner_join(results, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset1_Output/GLM_NBZI.txt", quote = F, sep="\t", row.names = F)
