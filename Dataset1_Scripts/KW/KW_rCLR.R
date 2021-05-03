### Perform differential abundance testing using Kruskal-Wallis rank-sum test on robust clr transformed abundances

date()

# Load required packages and source code
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
source("../../Support_Files/rclrMC.R")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Transform genus counts using robust clr with matrix completion
ps <- rclrMC(ps, ropt = NA, niter = 5, tol = 1e-5, seed=123)

# Loop Kruskal-Wallis test for all taxa
results <- data.frame()
for (i in 1:length(taxa_names(ps))){

  asv <- taxa_names(ps)[i]
 
  case_mclr <- mean(otu_table(subset_samples(ps, case_control == "Case"))[,i])
  cont_mclr <- mean(otu_table(subset_samples(ps, case_control == "Control"))[,i])
  mclrr <- case_mclr/cont_mclr
 
  pval <- kruskal.test(as.vector(otu_table(ps)[,i]) ~ sample_data(ps)$case_control)$p.value
 
  results <- rbind(results, data.frame(Representative_ASV = asv, 
                                       Case_Mean_CLR = case_mclr, 
                                       Control_Mean_CLR = cont_mclr, 
                                       Mean_CLR_Ratio = mclrr, 
                                       P = pval))
}

# Perform FDR correction for pvalues
results$FDR_BH <- p.adjust(results$P, method = 'BH')

# Add taxa designations to table
results <- inner_join(results, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset1_Output/KW_rCLR.txt", quote = F, sep="\t", row.names = F)
