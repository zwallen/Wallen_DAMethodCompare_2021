### Perform differential abundance testing using using Welch's t-test on relative abundances

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset2/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Log transform counts
ps <- transform_sample_counts(ps, function(x){log(x+1)})

# Transform log counts to relative abundance (total sum scaling)
ps <- transform_sample_counts(ps, function(x){x/sum(x)})

# Loop t-test for all taxa
results <- data.frame()
for (i in 1:length(taxa_names(ps))){

  asv <- taxa_names(ps)[i]
  
  case_mra <- mean(otu_table(subset_samples(ps, case_control == "Case"))[,i])
  cont_mra <- mean(otu_table(subset_samples(ps, case_control == "Control"))[,i])
  mrar <- case_mra/cont_mra
  
  pval <- t.test(as.vector(otu_table(ps)[,i]) ~ sample_data(ps)$case_control)$p.value
  
  results <- rbind(results, data.frame(Representative_ASV = asv, 
                                       Case_MRA = case_mra, 
                                       Control_MRA = cont_mra, 
                                       MRAR = mrar, 
                                       P = pval))
}
                                 
# Perform FDR correction for pvalues
results$FDR_BH <- p.adjust(results$P, method = 'BH')

# Add taxa designations to table
results <- inner_join(results, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset2_Output/t_test_TSS.txt", quote = F, sep="\t", row.names = F)
