### Perform differential abundance analysis using a generalized linear model on clr transformed abundances

date()

# Load required packages and source code
library(phyloseq); cat("Running phyloseq package version:", "\n"); packageVersion("phyloseq")
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")

# Read in phyloseq object
ps <- readRDS("../../PhyloseqObjects/Dataset1/phyloseq.rds")

# Collapse amplicon sequence variants to genus level
ps <- tax_glom(ps, taxrank = "Genus", NArm = F)

# Filter out taxa found in less than 10% of samples
ps <- filter_taxa(ps, function(x) sum(x > 0) >= (0.1*length(x)), TRUE)

# Transform genus counts using clr
ps <- transform_sample_counts(ps, function(x){log(x+1)-mean(log(x+1))})

# Dummy code case/control status
sample_data(ps)$case_control <- recode(sample_data(ps)$case_control, Case=1L, Control=0L)

# Apply GLM on clr transformed data for all taxa
results <- data.frame()
for (i in 1:length(taxa_names(ps))){

  asv <- taxa_names(ps)[i]
 
  mod <- glm(as.vector(otu_table(ps)[,i]) ~ sample_data(ps)$case_control)
  
  beta <- summary(mod)$coefficients[2,1]
  se <- summary(mod)$coefficients[2,2]
  pval <- summary(mod)$coefficients[2,4]

  results <- rbind(results, data.frame(Representative_ASV = asv, 
                                       Beta = beta,
				       SE = se,
				       P = pval))
}

# Perform FDR correction for pvalues
results$FDR_BH <- p.adjust(results$P, method = 'BH')

# Add taxa designations to table
results <- inner_join(results, rownames_to_column(data.frame(tax_table(ps)), "Representative_ASV"), by = "Representative_ASV")

write.table(results, "../../Script_Output/Dataset1_Output/GLM_CLR_10%.txt", quote = F, sep="\t", row.names = F)
