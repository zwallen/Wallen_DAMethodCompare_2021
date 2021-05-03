### Calculate pairwise concordances when taking replication across datasets into account

date()

# Load required packages
library(dplyr); cat("Running dplyr package version:", "\n"); packageVersion("dplyr")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")

##### Dataset 1 #####

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
#samseq <- read.table("../Script_Output/Dataset1_Output/SAMseq.txt", header=T, sep="\t")
lefse <- read.table("../Script_Output/Dataset1_Output/LEfSe.txt", header=T, sep="\t")

# Create input for calculating concordances (correlations)
input.1 <- data.frame(row.names = paste(kw_tss[order(kw_tss$Representative_ASV),]$Kingdom, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Phylum, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Class, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Order, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Family, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Genus, sep=";"))
input.1$`KW TSS`[kw_tss[order(kw_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`KW CLR`[kw_clr[order(kw_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`KW rCLR`[kw_rclr[order(kw_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`t-test TSS`[t_tss[order(t_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`t-test CLR`[t_clr[order(t_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`t-test rCLR`[t_rclr[order(t_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`GLM TSS`[glm_tss[order(glm_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`GLM CLR`[glm_clr[order(glm_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`GLM rCLR`[glm_rclr[order(glm_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`ANCOM`[ancom[order(ancom$Representative_ASV),]$detected_0.8 == "TRUE"] <- 1
input.1$`ANCOM-BC`[ancombc[order(ancombc$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`ALDEx2 t-test`[aldex2[order(aldex2$Representative_ASV),]$we.eBH < 0.05] <- 1
input.1$`ALDEx2 Wilcoxon`[aldex2[order(aldex2$Representative_ASV),]$wi.eBH < 0.05] <- 1
input.1$`baySeq`[bayseq[order(bayseq$Representative_ASV),]$FDR.DE < 0.05] <- 1
input.1$`DESeq2`[deseq2[order(deseq2$Representative_ASV),]$padj < 0.05] <- 1
input.1$`edgeR RLE`[edger_rle[order(edger_rle$Representative_ASV),]$FDR < 0.05] <- 1
input.1$`edgeR TMM`[edger_tmm[order(edger_tmm$Representative_ASV),]$FDR < 0.05] <- 1
input.1$`fitFeatureModel`[fitfeat[order(fitfeat$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.1$`fitZIG`[fitzig[order(fitzig$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.1$`GLM NBZI`[glm_nbzi[order(glm_nbzi$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.1$`limma-voom`[limma_voom[order(limma_voom$Representative_ASV),]$adj.P.Val < 0.05] <- 1
#input.1$`SAMseq`[samseq[order(samseq$Representative_ASV),]$q.value... < 5] <- 1
input.1$`LEfSe`[lefse[order(lefse$Representative_ASV),]$FDR_BH != "-"] <- 1
input.1[is.na(input.1)] <- 0
input.1$MRAR <- kw_tss[order(kw_tss$Representative_ASV),]$MRAR

##### Dataset 2 #####

# Read in results from methods
kw_tss <- read.table("../Script_Output/Dataset2_Output/KW_TSS.txt", header=T, sep="\t")
kw_clr <- read.table("../Script_Output/Dataset2_Output/KW_CLR.txt", header=T, sep="\t")
kw_rclr <- read.table("../Script_Output/Dataset2_Output/KW_rCLR.txt", header=T, sep="\t")
t_tss <- read.table("../Script_Output/Dataset2_Output/t_test_TSS.txt", header=T, sep="\t")
t_clr <- read.table("../Script_Output/Dataset2_Output/t_test_CLR.txt", header=T, sep="\t")
t_rclr <- read.table("../Script_Output/Dataset2_Output/t_test_rCLR.txt", header=T, sep="\t")
glm_tss <- read.table("../Script_Output/Dataset2_Output/GLM_TSS.txt", header=T, sep="\t")
glm_clr <- read.table("../Script_Output/Dataset2_Output/GLM_CLR.txt", header=T, sep="\t")
glm_rclr <- read.table("../Script_Output/Dataset2_Output/GLM_rCLR.txt", header=T, sep="\t")
ancom <- read.table("../Script_Output/Dataset2_Output/ANCOM.txt", header=T, sep="\t")
ancombc <- read.table("../Script_Output/Dataset2_Output/ANCOM-BC.txt", header=T, sep="\t")
aldex2 <- read.table("../Script_Output/Dataset2_Output/ALDEx2.txt", header=T, sep="\t")
bayseq <- read.table("../Script_Output/Dataset2_Output/baySeq.txt", header=T, sep="\t")
deseq2 <- read.table("../Script_Output/Dataset2_Output/DESeq2.txt", header=T, sep="\t")
edger_rle <- read.table("../Script_Output/Dataset2_Output/edgeR_RLE.txt", header=T, sep="\t")
edger_tmm <- read.table("../Script_Output/Dataset2_Output/edgeR_TMM.txt", header=T, sep="\t")
fitfeat <- read.table("../Script_Output/Dataset2_Output/fitFeatureModel.txt", header=T, sep="\t")
fitzig <- read.table("../Script_Output/Dataset2_Output/fitZIG.txt", header=T, sep="\t")
glm_nbzi <- read.table("../Script_Output/Dataset2_Output/GLM_NBZI.txt", header=T, sep="\t")
limma_voom <- read.table("../Script_Output/Dataset2_Output/limma_voom.txt", header=T, sep="\t")
#samseq <- read.table("../Script_Output/Dataset2_Output/SAMseq.txt", header=T, sep="\t")
lefse <- read.table("../Script_Output/Dataset2_Output/LEfSe.txt", header=T, sep="\t")

# Create input for calculating concordances (correlations)
input.2 <- data.frame(row.names = paste(kw_tss[order(kw_tss$Representative_ASV),]$Kingdom, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Phylum, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Class, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Order, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Family, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Genus, sep=";"))
input.2$`KW TSS`[kw_tss[order(kw_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`KW CLR`[kw_clr[order(kw_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`KW rCLR`[kw_rclr[order(kw_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`t-test TSS`[t_tss[order(t_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`t-test CLR`[t_clr[order(t_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`t-test rCLR`[t_rclr[order(t_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`GLM TSS`[glm_tss[order(glm_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`GLM CLR`[glm_clr[order(glm_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`GLM rCLR`[glm_rclr[order(glm_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`ANCOM`[ancom[order(ancom$Representative_ASV),]$detected_0.8 == "TRUE"] <- 1
input.2$`ANCOM-BC`[ancombc[order(ancombc$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`ALDEx2 t-test`[aldex2[order(aldex2$Representative_ASV),]$we.eBH < 0.05] <- 1
input.2$`ALDEx2 Wilcoxon`[aldex2[order(aldex2$Representative_ASV),]$wi.eBH < 0.05] <- 1
input.2$`baySeq`[bayseq[order(bayseq$Representative_ASV),]$FDR.DE < 0.05] <- 1
input.2$`DESeq2`[deseq2[order(deseq2$Representative_ASV),]$padj < 0.05] <- 1
input.2$`edgeR RLE`[edger_rle[order(edger_rle$Representative_ASV),]$FDR < 0.05] <- 1
input.2$`edgeR TMM`[edger_tmm[order(edger_tmm$Representative_ASV),]$FDR < 0.05] <- 1
input.2$`fitFeatureModel`[fitfeat[order(fitfeat$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.2$`fitZIG`[fitzig[order(fitzig$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.2$`GLM NBZI`[glm_nbzi[order(glm_nbzi$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.2$`limma-voom`[limma_voom[order(limma_voom$Representative_ASV),]$adj.P.Val < 0.05] <- 1
#input.2$`SAMseq`[samseq[order(samseq$Representative_ASV),]$q.value... < 5] <- 1
input.2$`LEfSe`[lefse[order(lefse$Representative_ASV),]$FDR_BH != "-"] <- 1
input.2[is.na(input.2)] <- 0
input.2$MRAR <- kw_tss[order(kw_tss$Representative_ASV),]$MRAR

##### Concordance calculation #####

# Filter for genera common between datasets
input.1 <- input.1[rownames(input.1) %in% rownames(input.2),]
input.2 <- input.2[rownames(input.2) %in% rownames(input.1),]

# Add matrices together to see who replicated based on significance and who did not
input <- input.1[order(rownames(input.1)),1:(ncol(input.1)-1)] + input.2[order(rownames(input.2)),1:(ncol(input.2)-1)]

# Differentiate replicated associations with same effect directions vs opposite effect directions
input[input == 2 & input.1[order(rownames(input.1)),]$MRAR > 1 & input.2[order(rownames(input.2)),]$MRAR > 1] <- 3
input[input == 2 & input.1[order(rownames(input.1)),]$MRAR < 1 & input.2[order(rownames(input.2)),]$MRAR < 1] <- 3

# Convert back to binary matrix
input[input != 3] <- 0
input[input == 3] <- 1

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
write.table(rownames_to_column(data.frame(conc.df[ncol(conc.df):1,], check.names=F), "Methods"), "../Script_Output/Joint_Analyses_Output/Pairwise_concordances.txt", row.names=F, quote=F, sep="\t")
