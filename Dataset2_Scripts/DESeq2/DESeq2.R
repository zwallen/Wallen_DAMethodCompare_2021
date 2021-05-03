### Perform differential abundance analysis using DESeq2

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(DESeq2); cat("Running DESeq2 package version:", "\n"); packageVersion("DESeq2")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset2/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Convert phyloseq object to DESeqDataSet object specifying design formula with only case/control status as predictor
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)
dds <- phyloseq_to_deseq2(ps, ~ case_control)

# Calculate size factors using modified geometric means that can handle zeros and estimate normalization factors
dds <- estimateSizeFactors(dds, type="poscounts")

# Perform differential abundance analysis
dds <- DESeq(dds)

# Get results and add taxonomy information
result.table <- results(dds, cooksCutoff = FALSE)
result.table <- inner_join(data.frame(Representative_ASV=rownames(result.table), result.table), rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(result.table, "../../Script_Output/Dataset2_Output/DESeq2.txt", row.names = F, quote = F, sep = '\t') 
