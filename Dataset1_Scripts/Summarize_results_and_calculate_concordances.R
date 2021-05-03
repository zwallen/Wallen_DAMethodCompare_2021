### Create summary of results and calculate pairwise concordances

date()

# Load required packages
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")

# Read in results from methods
kw_tss <- read.table("../Script_Output/Dataset1_Output/KW_TSS.txt", header=T, sep="\t")
kw_clr <- read.table("../Script_Output/Dataset1_Output/KW_CLR.txt", header=T, sep="\t")
kw_rclr <- read.table("../Script_Output/Dataset1_Output/KW_rCLR.txt", header=T, sep="\t")
t_tss <- read.table("../Script_Output/Dataset1_Output/t_test_TSS.txt", header=T, sep="\t")
t_clr <- read.table("../Script_Output/Dataset1_Output/t_test_CLR.txt", header=T, sep="\t")
t_rclr <- read.table("../Script_Output/Dataset1_Output/t_test_rCLR.txt", header=T, sep="\t")
glm_tss <- read.table("../Script_Output/Dataset1_Output/GLM_TSS.txt", header=T, sep="\t")
glm_clr <- read.table("../Script_Output/Dataset1_Output/GLM_CLR.txt", header=T, sep="\t")
glm_rclr <- read.table("../Script_Output/Dataset1_Output/GLM_rCLR.txt", header=T, sep="\t")
ancom <- read.table("../Script_Output/Dataset1_Output/ANCOM.txt", header=T, sep="\t")
ancombc <- read.table("../Script_Output/Dataset1_Output/ANCOM-BC.txt", header=T, sep="\t")
aldex2 <- read.table("../Script_Output/Dataset1_Output/ALDEx2.txt", header=T, sep="\t")
bayseq <- read.table("../Script_Output/Dataset1_Output/baySeq.txt", header=T, sep="\t")
deseq2 <- read.table("../Script_Output/Dataset1_Output/DESeq2.txt", header=T, sep="\t")
edger_rle <- read.table("../Script_Output/Dataset1_Output/edgeR_RLE.txt", header=T, sep="\t")
edger_tmm <- read.table("../Script_Output/Dataset1_Output/edgeR_TMM.txt", header=T, sep="\t")
fitfeat <- read.table("../Script_Output/Dataset1_Output/fitFeatureModel.txt", header=T, sep="\t")
fitzig <- read.table("../Script_Output/Dataset1_Output/fitZIG.txt", header=T, sep="\t")
glm_nbzi <- read.table("../Script_Output/Dataset1_Output/GLM_NBZI.txt", header=T, sep="\t")
limma_voom <- read.table("../Script_Output/Dataset1_Output/limma_voom.txt", header=T, sep="\t")
samseq <- read.table("../Script_Output/Dataset1_Output/SAMseq.txt", header=T, sep="\t")
lefse <- read.table("../Script_Output/Dataset1_Output/LEfSe.txt", header=T, sep="\t")

# Combine FDR qvalues into one table
results_all <- data.frame(`KW TSS` = kw_tss[order(kw_tss$Representative_ASV),]$FDR_BH,
                          `KW CLR` = kw_clr[order(kw_clr$Representative_ASV),]$FDR_BH,
                          `KW rCLR` = kw_rclr[order(kw_rclr$Representative_ASV),]$FDR_BH,
                          `t-test TSS` = t_tss[order(t_tss$Representative_ASV),]$FDR_BH,
                          `t-test CLR` = t_clr[order(t_clr$Representative_ASV),]$FDR_BH,
                          `t-test rCLR` = t_rclr[order(t_rclr$Representative_ASV),]$FDR_BH,
                          `GLM TSS` = glm_tss[order(glm_tss$Representative_ASV),]$FDR_BH,
                          `GLM CLR` = glm_clr[order(glm_clr$Representative_ASV),]$FDR_BH,
                          `GLM rCLR` = glm_rclr[order(glm_rclr$Representative_ASV),]$FDR_BH,
                          `ANCOM` = abs(as.integer(ancom[order(ancom$Representative_ASV),]$detected_0.8)-1),
                          `ANCOM-BC` = ancombc[order(ancombc$Representative_ASV),]$FDR_BH,
                          `ALDEx2 t-test` = aldex2[order(aldex2$Representative_ASV),]$we.eBH,
                          `ALDEx2 Wilcoxon` = aldex2[order(aldex2$Representative_ASV),]$wi.eBH,
                          `baySeq` = bayseq[order(bayseq$Representative_ASV),]$FDR.DE,
                          `DESeq2` = deseq2[order(deseq2$Representative_ASV),]$padj,
                          `edgeR RLE` = edger_rle[order(edger_rle$Representative_ASV),]$FDR,
                          `edgeR TMM` = edger_tmm[order(edger_tmm$Representative_ASV),]$FDR,
                          `fitFeatureModel` = fitfeat[order(fitfeat$Representative_ASV),]$adjPvalues,
                          `fitZIG` = fitzig[order(fitzig$Representative_ASV),]$adjPvalues,
                          `GLM NBZI` = glm_nbzi[order(glm_nbzi$Representative_ASV),]$FDR_BH,
                          `limma-voom` = limma_voom[order(limma_voom$Representative_ASV),]$adj.P.Val,
                          `SAMseq` = (samseq[order(samseq$Representative_ASV),]$q.value...)/100,
                          `LEfSe` = as.numeric(lefse[order(lefse$Representative_ASV),]$FDR_BH), check.names=F)
results_all[is.na(results_all)] <- 1

# Calculate the number of methods that detected each taxon as significantly different between cases and controls
results_all$`N methods detected by` <- rowSums(results_all < 0.05)

# Calculate the number and proportion of taxa each method detected as significantly different between cases and controls
results_all <- rbind(results_all, colSums(results_all < 0.05), colSums(results_all < 0.05)/nrow(results_all))

# Add taxonomic information and rownames for bottom two summary rows, then write out
results_all <- data.frame(`Kingdom` = c(kw_tss[order(kw_tss$Representative_ASV),]$Kingdom, " ", " "),
                          `Phylum` = c(kw_tss[order(kw_tss$Representative_ASV),]$Phylum, " ", " "),
                          `Class` = c(kw_tss[order(kw_tss$Representative_ASV),]$Class, " ", " "),
                          `Order` = c(kw_tss[order(kw_tss$Representative_ASV),]$Order, " ", " "),
                          `Family` = c(kw_tss[order(kw_tss$Representative_ASV),]$Family, " ", " "),
                          `Genus` = c(kw_tss[order(kw_tss$Representative_ASV),]$Genus, "N significant", "Proportion significant"),
                          `Case MRA` = c(kw_tss[order(kw_tss$Representative_ASV),]$Case_MRA, " ", " "),
                          `Control MRA` = c(kw_tss[order(kw_tss$Representative_ASV),]$Control_MRA, " ", " "),
                          `MRAR` = c(kw_tss[order(kw_tss$Representative_ASV),]$MRAR, " ", " "),
                          results_all, check.names=F)
results_all[c(nrow(results_all)-1,nrow(results_all)), ncol(results_all)] <- " "
write.table(results_all, "../Script_Output/Dataset1_Output/Results_summary_table.txt", row.names=F, quote=F, sep="\t")

# Create input for calculating concordances (correlations)
input <- data.frame(row.names = 1:dim(kw_tss)[1])
input$`KW TSS`[kw_tss[order(kw_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`KW CLR`[kw_clr[order(kw_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`KW rCLR`[kw_rclr[order(kw_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`t-test TSS`[t_tss[order(t_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`t-test CLR`[t_clr[order(t_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`t-test rCLR`[t_rclr[order(t_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`GLM TSS`[glm_tss[order(glm_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`GLM CLR`[glm_clr[order(glm_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`GLM rCLR`[glm_rclr[order(glm_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`ANCOM`[ancom[order(ancom$Representative_ASV),]$detected_0.8 == "TRUE"] <- 1
input$`ANCOM-BC`[ancombc[order(ancombc$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`ALDEx2 t-test`[aldex2[order(aldex2$Representative_ASV),]$we.eBH < 0.05] <- 1
input$`ALDEx2 Wilcoxon`[aldex2[order(aldex2$Representative_ASV),]$wi.eBH < 0.05] <- 1
input$`baySeq`[bayseq[order(bayseq$Representative_ASV),]$FDR.DE < 0.05] <- 1
input$`DESeq2`[deseq2[order(deseq2$Representative_ASV),]$padj < 0.05] <- 1
input$`edgeR RLE`[edger_rle[order(edger_rle$Representative_ASV),]$FDR < 0.05] <- 1
input$`edgeR TMM`[edger_tmm[order(edger_tmm$Representative_ASV),]$FDR < 0.05] <- 1
input$`fitFeatureModel`[fitfeat[order(fitfeat$Representative_ASV),]$adjPvalues < 0.05] <- 1
input$`fitZIG`[fitzig[order(fitzig$Representative_ASV),]$adjPvalues < 0.05] <- 1
input$`GLM NBZI`[glm_nbzi[order(glm_nbzi$Representative_ASV),]$FDR_BH < 0.05] <- 1
input$`limma-voom`[limma_voom[order(limma_voom$Representative_ASV),]$adj.P.Val < 0.05] <- 1
input$`SAMseq`[samseq[order(samseq$Representative_ASV),]$q.value... < 5] <- 1
input$`LEfSe`[lefse[order(lefse$Representative_ASV),]$FDR_BH != "-"] <- 1
input[is.na(input)] <- 0

# Calculate pairwise concordances for method results
conc.df <- apply(input, 2, 
                 function(x){apply(input, 2, 
                                   function(y){sum(x + y == 2)/sum(x + y > 0)})
                             }
                 )

# Order methods by mean concordance
means <- sort(apply(conc.df, 2, mean))
conc.df <- conc.df[match(names(means), rownames(conc.df)), match(names(means), colnames(conc.df))]

# Write table of concordance values out
write.table(rownames_to_column(data.frame(conc.df[ncol(conc.df):1,], check.names=F), "Methods"), "../Script_Output/Dataset1_Output/Pairwise_concordances.txt", row.names=F, quote=F, sep="\t")
