### Perform differential abundance testing using zero-inflated lognormal model

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(metagenomeSeq); cat("Running metagenomeSeq package version:", "\n"); packageVersion("metagenomeSeq")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset2/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Convert phyloseq object to MRexperiment object
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)
mre <- phyloseq_to_metagenomeSeq(ps)

# Perform Cummulative Sum Scaling transformation on data
mre <- cumNorm(mre)

# Perform differential abundance analysis
mod <- with(pData(mre), model.matrix(~ case_control))
fit.results <- fitFeatureModel(mre, mod = mod)

# FDR correct Pvalues and extract results adding taxonomy information to result table
result.table <- MRfulltable(fit.results, number = nrow(assayData(mre)$counts), adjustMethod = "BH")
result.table <- inner_join(rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"),
                           rownames_to_column(result.table, "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(result.table, "../../Script_Output/Dataset2_Output/fitFeatureModel_10%.txt", row.names = F, quote = F, sep = '\t')
