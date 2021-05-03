### Perform differential abundance analysis using ANCOM-BC

date()

# Load required packages and source code
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(ANCOMBC); cat("Running ANCOMBC package version:", "\n"); packageVersion("ANCOMBC")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Dummy code case/control status
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)

# Run ANCOM-BC
ancom.out <- ancombc(phyloseq = ps,
                     formula = "case_control",
                     p_adj_method = "BH",
                     zero_cut = 1,
                     lib_cut = 0,
                     group = NULL,
                     struc_zero = FALSE,
                     neg_lb = FALSE,
                     tol = 1e-05,
                     max_iter = 100,
                     conserve = FALSE,
                     alpha = 0.05,
                     global = FALSE)

# Extract results and add taxa designations to table
ancom.res <- data.frame(Representative_ASV=rownames(ancom.out$res$beta), Beta=ancom.out$res$beta[,1], 
                        SE=ancom.out$res$se[,1], P=ancom.out$res$p_val[,1], FDR_BH=ancom.out$res$q_val[,1])
results <- inner_join(ancom.res, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

# Write out results
write.table(results, "../../Script_Output/Dataset1_Output/ANCOM-BC_10%.txt", quote = F, sep="\t", row.names = F)
