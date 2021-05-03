### Perform differential abundance testing using using limma linear modeling with voom transformation and weights

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(edgeR); cat("Running edgeR package version:", "\n"); packageVersion("edgeR")
library(limma); cat("Running limma package version:", "\n"); packageVersion("limma")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Convert phyloseq object to DGEList object specifying groups to be tested
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)
dge <- DGEList(counts=t(as(otu_table(ps), "matrix")), 
		group=sample_data(ps)$case_control, 
		genes=tax_table(ps))
		
# Calculate tagwise dispersions
dge <- calcNormFactors(dge)

# Create model matrix
mod <- with(data.frame(sample_data(ps)), model.matrix(~ case_control))

# Use function voom to transform to log2-CPM and calculate associated weights to ready data for linear modeling with limma
dge.voom <- voom(dge, design=mod, save.plot=T)

# Perform linear modeling and calculate empirical Bayes statistics for differential expression
fit <- lmFit(dge.voom, design=mod)
fit <- eBayes(fit)

# FDR correct Pvalues and extract results
result.table <- topTable(fit, coef="case_control", n = Inf, sort.by = "p", adjust.method="BH")
result.table <- data.frame(Representative_ASV=rownames(result.table), result.table)

# Write out results
write.table(result.table, "../../Script_Output/Dataset1_Output/limma_voom_10%.txt", row.names = F, quote = F, sep = '\t') 
