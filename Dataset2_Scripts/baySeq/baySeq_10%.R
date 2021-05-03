### Perform differential abundance testing using baySeq

date()

# Load required packages
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(baySeq); cat("Running baySeq package version:", "\n"); packageVersion("baySeq")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset2/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Calculate total sequence number to use as size factor for baySeq
sample_data(ps)$total_seq_count <- rowSums(otu_table(ps))

# Combine data into countData object and supply library size
replicates <- sample_data(ps)$case_control
groups <- list(NDE = rep(0, dim(otu_table(ps))[1]), DE = recode(sample_data(ps)$case_control, Case=1L, Control=0L))
count.data <- t(data.frame(otu_table(ps)))
cd <- new("countData", data=count.data, groups=groups, replicates=replicates)
libsizes(cd) <- sample_data(ps)$total_seq_count

# Estimate prior parameters using default negative binomial distribution
cd <- getPriors.NB(cd, cl = NULL)

# Find posterior likelihood of each count or paired count belonging to one of the defined models
cd <- getLikelihoods(cd, cl = NULL)

# Get results
results <- topCounts(cd, group = "DE", number = Inf)

# Add taxa designations to table
results <- inner_join(rownames_to_column(results[,(dim(cd@data)[2]+2):dim(results)[2]], "Representative_ASV"),
		      rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"),
                      by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset2_Output/baySeq_10%.txt", quote = F, sep="\t", row.names = F)
