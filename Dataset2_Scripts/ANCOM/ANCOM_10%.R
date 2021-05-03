### Perform differential abundance analysis using ANCOM (updated R code)

date()

# Load required packages and source code
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
source("../../Support_Files/ANCOM_updated_code.R")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset2/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Extract taxa count data from phyloseq object
taxa.df <- data.frame(otu_table(ps))
taxa.df <- rownames_to_column(taxa.df, "Sample.ID")

# Extract sample metadata from phyloseq object
meta.df <- data.frame(sample_data(ps))
meta.df <- meta.df[,-1]
colnames(meta.df)[1] <- "Sample.ID"

# Run ANCOM
ancom.results <- ANCOM.main(OTUdat=taxa.df,
                           Vardat=meta.df,
                           adjusted=FALSE,
                           repeated=FALSE,
                           main.var="case_control", 
                           adj.formula=NULL,
                           repeat.var=NULL,
                           longitudinal=FALSE,
                           random.formula=NULL,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=1)

# Extract results and add taxa designations to table
ancom.Wtaxa <- ancom.results$W.taxa
colnames(ancom.Wtaxa)[1] <- "Representative_ASV"
results <- inner_join(ancom.Wtaxa, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset2_Output/ANCOM_10%.txt", quote = F, sep="\t", row.names = F)
