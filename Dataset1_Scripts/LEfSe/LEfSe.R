### Perform differential abundance testing using LEfSe

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Transform counts to relative abundance (total sum scaling)
ps <- transform_sample_counts(ps, function(x){x/sum(x)})

# Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps))
taxa.df <- rownames_to_column(taxa.df, "SampleID")

# Add column for case control designations
taxa.df <- cbind(data.frame(case_control=sample_data(ps)$case_control), taxa.df)

# Format for LEfSe input
lefse.df <- rownames_to_column(as.data.frame(t(taxa.df)))

# Write out LEfSe input data
write.table(lefse.df, "../../Script_Output/Dataset1_Output/LEfSe_input.txt", quote = F, sep="\t", row.names = F, col.names = F)

### Run LEfSe workflow ###
# Requires for LEfSe scripts and their dependencies to be in your PATH, which I found is easiest by making a conda environment
# for LEfSe with Python 2 and downloading through conda

# Run python script for formatting
system("lefse-format_input.py ../../Script_Output/Dataset1_Output/LEfSe_input.txt ../../Script_Output/Dataset1_Output/LEfSe_input.in -c 1 -s -1 -u 2 -o 1000000")

# Run LEfSe using default parameters
system("run_lefse.py ../../Script_Output/Dataset1_Output/LEfSe_input.in ../../Script_Output/Dataset1_Output/LEfSe.txt")

# Run LEfSe using parameters that will output all LDA scores and pvalues
system("run_lefse.py ../../Script_Output/Dataset1_Output/LEfSe_input.in ../../Script_Output/Dataset1_Output/LEfSe_full.txt -a 1.0 -w 1.0 -l 0.0")

# Clean up input files
system("rm ../../Script_Output/Dataset1_Output/LEfSe_input.txt")
system("rm ../../Script_Output/Dataset1_Output/LEfSe_input.in")

##########################

# Read in results
lefse.results <- read.table("../../Script_Output/Dataset1_Output/LEfSe.txt", sep="\t")
lefse.results.full <- read.table("../../Script_Output/Dataset1_Output/LEfSe_full.txt")
colnames(lefse.results) <- c("Representative_ASV", "LogOfHighestClassAvg", "HighestClass", "LDAScore", "P")

# Add FDR q-values in LEfSe results and further filter for those that are FDR < 0.05
lefse.results$FDR_BH <- p.adjust(lefse.results.full[,5], method="BH")
lefse.results[lefse.results$FDR_BH > 0.05,c(5,6)] <- "-"
lefse.results[lefse.results$P == "-",6] <- "-"
lefse.results[is.na(lefse.results$LDAScore),6] <- "-"
lefse.results[lefse.results[,6] == "-",3] <- ""
lefse.results[lefse.results[,6] == "-",4] <- NA

# Add taxa designations to table
lefse.results <- inner_join(lefse.results, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results now with corrected pvalues incorporated
write.table(lefse.results, "../../Script_Output/Dataset1_Output/LEfSe.txt", quote = F, sep="\t", row.names = F)
system("rm ../../Script_Output/Dataset1_Output/LEfSe_full.txt")
