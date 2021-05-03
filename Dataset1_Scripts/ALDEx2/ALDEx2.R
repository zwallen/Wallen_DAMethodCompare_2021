### Perform differential abundance testing using ALDEx2 Welch's t-test and Wilcoxon test

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(ALDEx2); cat("Running ALDEx2 package version:", "\n"); packageVersion("ALDEx2")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Set the comparison groups
conds <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)

# Perform ALDEx2 analysis
set.seed(1234)
x <- aldex(t(data.frame(otu_table(ps))), 
	   conds, 
	   mc.samples=1000,
	   test="t",
	   effect=TRUE,
	   include.sample.summary=F,
	   denom="all",
	   verbose=TRUE)

# Add taxa designations to table
x.all <- inner_join(rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"),
		    rownames_to_column(x, "Representative_ASV"),
		    by = "Representative_ASV")

# Write out results
write.table(x.all, "../../Script_Output/Dataset1_Output/ALDEx2.txt", quote = F, sep="\t", row.names = F)
