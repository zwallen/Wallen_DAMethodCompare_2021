### Perform differential abundance testing using SAMseq

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(samr); cat("Running samr package version:", "\n"); packageVersion("samr")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Recode case/control variable with 1 and 2 (required by SAMseq)
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=2L, Control=1L)

# Normalize data using SAM method and round to whole integer (required by SAMseq)
norm.data <- round(samr.norm.data(t(otu_table(ps))), 0)

# Run SAMseq analysis
set.seed(1234)
sam.results <- SAMseq(norm.data, sample_data(ps)$case_control, resp.type = "Two class unpaired", geneid = taxa_names(ps), fdr.output = 1)

# Create results table
results <- data.frame(Representative_ASV = rownames(sam.results$samr.obj$x), FoldChange = sam.results$samr.obj$foldchange) 
sig.gene.table <- rbind(sam.results$siggenes.table$genes.up, sam.results$siggenes.table$genes.lo)[,-1]
colnames(sig.gene.table)[1] <- "Representative_ASV"
results <- full_join(data.frame(results), data.frame(sig.gene.table), by = "Representative_ASV")

# Add taxa designations to results
results <- inner_join(results, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset1_Output/SAMseq_10%.txt", quote = F, sep="\t", row.names = F)
