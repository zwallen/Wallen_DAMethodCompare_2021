### Perform differential abundance testing using edgeR exactTest with TMM calculated scaling factors

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(edgeR); cat("Running edgeR package version:", "\n"); packageVersion("edgeR")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Convert phyloseq object to DGEList object specifying groups to be tested
dge <- DGEList(counts=t(as(otu_table(ps), "matrix")), 
		group=recode(sample_data(ps)$case_control, Case=1L, Control=0L), 
		genes=tax_table(ps))

# Calculate normalization factors and tagwise dispersions
dge <- calcNormFactors(dge)
dge <- estimateTagwiseDisp(estimateCommonDisp(dge))

# Perform differential abundance analysis
et.results <- exactTest(dge)

# FDR correct Pvalues and extract results
result.table <- topTags(et.results, n = nrow(dge$table), adjust.method = "BH", sort.by = "PValue")
result.table <- data.frame(Representative_ASV=rownames(result.table), result.table)

# Write out results
write.table(result.table, "../../Script_Output/Dataset1_Output/edgeR_TMM_10%.txt", row.names = F, quote = F, sep = '\t') 
