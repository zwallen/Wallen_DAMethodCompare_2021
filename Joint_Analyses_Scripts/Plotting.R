### Create plots used for figures in the manuscript

date()

# Load required packages and source code
library(reshape2); cat("Running reshape2 package version:", "\n"); packageVersion("reshape2")
library(tibble); cat("Running tibble package version:", "\n"); packageVersion("tibble")
library(cluster); cat("Running cluster package version:", "\n"); packageVersion("cluster")
library(ggplot2); cat("Running ggplot2 package version:", "\n"); packageVersion("ggplot2")
library(ggh4x); cat("Running ggh4x package version:", "\n"); packageVersion("ggh4x")
library(ggpubr); cat("Running ggpubr package version:", "\n"); packageVersion("ggpubr")
library(grid); cat("Running grid package version:", "\n"); packageVersion("grid")
library(scales); cat("Running scales package version:", "\n"); packageVersion("scales")
source("../Support_Files/heatmap.3.R")

###########################################################################
##### Create plot summarizing pairwise concordances and proportion of #####
##### DA signatures detected by methods within datasets and for       #####
##### replicated signatures                                           #####
###########################################################################

##### Prepping proportion of detected differentially abundant taxa data #####

##### Dataset 1 unfiltered #####

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

# Create binary taxa by method matrix
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
input.1$`SAMseq`[samseq[order(samseq$Representative_ASV),]$q.value... < 5] <- 1
input.1$`LEfSe`[lefse[order(lefse$Representative_ASV),]$FDR_BH != "-"] <- 1
input.1[is.na(input.1)] <- 0
input.1$MRA <- (kw_tss[order(kw_tss$Representative_ASV),]$Case_MRA +
                kw_tss[order(kw_tss$Representative_ASV),]$Control_MRA)/2
input.1$MRAR <- kw_tss[order(kw_tss$Representative_ASV),]$MRAR

##### Dataset 2 unfiltered #####

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

# Create binary taxa by method matrix
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
input.2$MRA <- (kw_tss[order(kw_tss$Representative_ASV),]$Case_MRA +
                kw_tss[order(kw_tss$Representative_ASV),]$Control_MRA)/2
input.2$MRAR <- kw_tss[order(kw_tss$Representative_ASV),]$MRAR

##### Replicated signatures unfiltered #####

# Filter for genera common between datasets
input.1.sub <- input.1[rownames(input.1) %in% rownames(input.2),-(ncol(input.1)-3)]
input.2.sub <- input.2[rownames(input.2) %in% rownames(input.1),]

# Add matrices together to see who replicated based on significance and who did not
input.3 <- input.1.sub[order(rownames(input.1.sub)),1:(ncol(input.1.sub)-2)] + input.2.sub[order(rownames(input.2.sub)),1:(ncol(input.2.sub)-2)]

# Differentiate replicated associations with same effect directions vs opposite effect directions
input.3[input.3 == 2 & input.1.sub[order(rownames(input.1.sub)),]$MRAR > 1 & input.2.sub[order(rownames(input.2.sub)),]$MRAR > 1] <- 3
input.3[input.3 == 2 & input.1.sub[order(rownames(input.1.sub)),]$MRAR < 1 & input.2.sub[order(rownames(input.2.sub)),]$MRAR < 1] <- 3

# Convert back to binary matrix
input.3[input.3 != 3] <- 0
input.3[input.3 == 3] <- 1

# Add MRAs and MRARs from both datasets
input.3 <- cbind(input.3, data.frame(MRAR_1=input.1.sub[order(rownames(input.1.sub)),]$MRAR,
                                     MRAR_2=input.2.sub[order(rownames(input.2.sub)),]$MRAR,
                                     MRA_1=input.1.sub[order(rownames(input.1.sub)),]$MRA,
                                     MRA_2=input.2.sub[order(rownames(input.2.sub)),]$MRA))

# Print N replicated
rep1 <- data.frame(Taxa=names(rowSums(input.3[,1:(ncol(input.3)-4)])), 
                   `N methods replicated`=rowSums(input.3[,1:(ncol(input.3)-4)]), check.names=F)
write.table(rep1, "../Script_Output/Joint_Analyses_Output/N_replicated.txt", quote=F, row.names=F, sep="\t")

##### Dataset 1 filtered #####

# Read in results from methods
kw_tss <- read.table("../Script_Output/Dataset1_Output/KW_TSS_10%.txt", header=T, sep="\t")
kw_clr <- read.table("../Script_Output/Dataset1_Output/KW_CLR_10%.txt", header=T, sep="\t")
kw_rclr <- read.table("../Script_Output/Dataset1_Output/KW_rCLR_10%.txt", header=T, sep="\t")
t_tss <- read.table("../Script_Output/Dataset1_Output/t_test_TSS_10%.txt", header=T, sep="\t")
t_clr <- read.table("../Script_Output/Dataset1_Output/t_test_CLR_10%.txt", header=T, sep="\t")
t_rclr <- read.table("../Script_Output/Dataset1_Output/t_test_rCLR_10%.txt", header=T, sep="\t")
glm_tss <- read.table("../Script_Output/Dataset1_Output/GLM_TSS_10%.txt", header=T, sep="\t")
glm_clr <- read.table("../Script_Output/Dataset1_Output/GLM_CLR_10%.txt", header=T, sep="\t")
glm_rclr <- read.table("../Script_Output/Dataset1_Output/GLM_rCLR_10%.txt", header=T, sep="\t")
ancom <- read.table("../Script_Output/Dataset1_Output/ANCOM_10%.txt", header=T, sep="\t")
ancombc <- read.table("../Script_Output/Dataset1_Output/ANCOM-BC_10%.txt", header=T, sep="\t")
aldex2 <- read.table("../Script_Output/Dataset1_Output/ALDEx2_10%.txt", header=T, sep="\t")
bayseq <- read.table("../Script_Output/Dataset1_Output/baySeq_10%.txt", header=T, sep="\t")
deseq2 <- read.table("../Script_Output/Dataset1_Output/DESeq2_10%.txt", header=T, sep="\t")
edger_rle <- read.table("../Script_Output/Dataset1_Output/edgeR_RLE_10%.txt", header=T, sep="\t")
edger_tmm <- read.table("../Script_Output/Dataset1_Output/edgeR_TMM_10%.txt", header=T, sep="\t")
fitfeat <- read.table("../Script_Output/Dataset1_Output/fitFeatureModel_10%.txt", header=T, sep="\t")
fitzig <- read.table("../Script_Output/Dataset1_Output/fitZIG_10%.txt", header=T, sep="\t")
glm_nbzi <- read.table("../Script_Output/Dataset1_Output/GLM_NBZI_10%.txt", header=T, sep="\t")
limma_voom <- read.table("../Script_Output/Dataset1_Output/limma_voom_10%.txt", header=T, sep="\t")
samseq <- read.table("../Script_Output/Dataset1_Output/SAMseq_10%.txt", header=T, sep="\t")
lefse <- read.table("../Script_Output/Dataset1_Output/LEfSe_10%.txt", header=T, sep="\t")

# Create binary taxa by method matrix
input.4 <- data.frame(row.names = paste(kw_tss[order(kw_tss$Representative_ASV),]$Kingdom, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Phylum, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Class, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Order, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Family, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Genus, sep=";"))
input.4$`KW TSS`[kw_tss[order(kw_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`KW CLR`[kw_clr[order(kw_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`KW rCLR`[kw_rclr[order(kw_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`t-test TSS`[t_tss[order(t_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`t-test CLR`[t_clr[order(t_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`t-test rCLR`[t_rclr[order(t_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`GLM TSS`[glm_tss[order(glm_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`GLM CLR`[glm_clr[order(glm_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`GLM rCLR`[glm_rclr[order(glm_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`ANCOM`[ancom[order(ancom$Representative_ASV),]$detected_0.8 == "TRUE"] <- 1
input.4$`ANCOM-BC`[ancombc[order(ancombc$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`ALDEx2 t-test`[aldex2[order(aldex2$Representative_ASV),]$we.eBH < 0.05] <- 1
input.4$`ALDEx2 Wilcoxon`[aldex2[order(aldex2$Representative_ASV),]$wi.eBH < 0.05] <- 1
input.4$`baySeq`[bayseq[order(bayseq$Representative_ASV),]$FDR.DE < 0.05] <- 1
input.4$`DESeq2`[deseq2[order(deseq2$Representative_ASV),]$padj < 0.05] <- 1
input.4$`edgeR RLE`[edger_rle[order(edger_rle$Representative_ASV),]$FDR < 0.05] <- 1
input.4$`edgeR TMM`[edger_tmm[order(edger_tmm$Representative_ASV),]$FDR < 0.05] <- 1
input.4$`fitFeatureModel`[fitfeat[order(fitfeat$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.4$`fitZIG`[fitzig[order(fitzig$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.4$`GLM NBZI`[glm_nbzi[order(glm_nbzi$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.4$`limma-voom`[limma_voom[order(limma_voom$Representative_ASV),]$adj.P.Val < 0.05] <- 1
input.4$`SAMseq`[samseq[order(samseq$Representative_ASV),]$q.value... < 5] <- 1
input.4$`LEfSe`[lefse[order(lefse$Representative_ASV),]$FDR_BH != "-"] <- 1
input.4[is.na(input.4)] <- 0
input.4$MRA <- (kw_tss[order(kw_tss$Representative_ASV),]$Case_MRA +
                kw_tss[order(kw_tss$Representative_ASV),]$Control_MRA)/2
input.4$MRAR <- kw_tss[order(kw_tss$Representative_ASV),]$MRAR

##### Dataset 2 filtered #####

# Read in results from methods
kw_tss <- read.table("../Script_Output/Dataset2_Output/KW_TSS_10%.txt", header=T, sep="\t")
kw_clr <- read.table("../Script_Output/Dataset2_Output/KW_CLR_10%.txt", header=T, sep="\t")
kw_rclr <- read.table("../Script_Output/Dataset2_Output/KW_rCLR_10%.txt", header=T, sep="\t")
t_tss <- read.table("../Script_Output/Dataset2_Output/t_test_TSS_10%.txt", header=T, sep="\t")
t_clr <- read.table("../Script_Output/Dataset2_Output/t_test_CLR_10%.txt", header=T, sep="\t")
t_rclr <- read.table("../Script_Output/Dataset2_Output/t_test_rCLR_10%.txt", header=T, sep="\t")
glm_tss <- read.table("../Script_Output/Dataset2_Output/GLM_TSS_10%.txt", header=T, sep="\t")
glm_clr <- read.table("../Script_Output/Dataset2_Output/GLM_CLR_10%.txt", header=T, sep="\t")
glm_rclr <- read.table("../Script_Output/Dataset2_Output/GLM_rCLR_10%.txt", header=T, sep="\t")
ancom <- read.table("../Script_Output/Dataset2_Output/ANCOM_10%.txt", header=T, sep="\t")
ancombc <- read.table("../Script_Output/Dataset2_Output/ANCOM-BC_10%.txt", header=T, sep="\t")
aldex2 <- read.table("../Script_Output/Dataset2_Output/ALDEx2_10%.txt", header=T, sep="\t")
bayseq <- read.table("../Script_Output/Dataset2_Output/baySeq_10%.txt", header=T, sep="\t")
deseq2 <- read.table("../Script_Output/Dataset2_Output/DESeq2_10%.txt", header=T, sep="\t")
edger_rle <- read.table("../Script_Output/Dataset2_Output/edgeR_RLE_10%.txt", header=T, sep="\t")
edger_tmm <- read.table("../Script_Output/Dataset2_Output/edgeR_TMM_10%.txt", header=T, sep="\t")
fitfeat <- read.table("../Script_Output/Dataset2_Output/fitFeatureModel_10%.txt", header=T, sep="\t")
fitzig <- read.table("../Script_Output/Dataset2_Output/fitZIG_10%.txt", header=T, sep="\t")
glm_nbzi <- read.table("../Script_Output/Dataset2_Output/GLM_NBZI_10%.txt", header=T, sep="\t")
limma_voom <- read.table("../Script_Output/Dataset2_Output/limma_voom_10%.txt", header=T, sep="\t")
samseq <- read.table("../Script_Output/Dataset2_Output/SAMseq_10%.txt", header=T, sep="\t")
lefse <- read.table("../Script_Output/Dataset2_Output/LEfSe_10%.txt", header=T, sep="\t")

# Create binary taxa by method matrix
input.5 <- data.frame(row.names = paste(kw_tss[order(kw_tss$Representative_ASV),]$Kingdom, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Phylum, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Class, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Order, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Family, 
                                        kw_tss[order(kw_tss$Representative_ASV),]$Genus, sep=";"))
input.5$`KW TSS`[kw_tss[order(kw_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`KW CLR`[kw_clr[order(kw_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`KW rCLR`[kw_rclr[order(kw_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`t-test TSS`[t_tss[order(t_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`t-test CLR`[t_clr[order(t_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`t-test rCLR`[t_rclr[order(t_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`GLM TSS`[glm_tss[order(glm_tss$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`GLM CLR`[glm_clr[order(glm_clr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`GLM rCLR`[glm_rclr[order(glm_rclr$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`ANCOM`[ancom[order(ancom$Representative_ASV),]$detected_0.8 == "TRUE"] <- 1
input.5$`ANCOM-BC`[ancombc[order(ancombc$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`ALDEx2 t-test`[aldex2[order(aldex2$Representative_ASV),]$we.eBH < 0.05] <- 1
input.5$`ALDEx2 Wilcoxon`[aldex2[order(aldex2$Representative_ASV),]$wi.eBH < 0.05] <- 1
input.5$`baySeq`[bayseq[order(bayseq$Representative_ASV),]$FDR.DE < 0.05] <- 1
input.5$`DESeq2`[deseq2[order(deseq2$Representative_ASV),]$padj < 0.05] <- 1
input.5$`edgeR RLE`[edger_rle[order(edger_rle$Representative_ASV),]$FDR < 0.05] <- 1
input.5$`edgeR TMM`[edger_tmm[order(edger_tmm$Representative_ASV),]$FDR < 0.05] <- 1
input.5$`fitFeatureModel`[fitfeat[order(fitfeat$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.5$`fitZIG`[fitzig[order(fitzig$Representative_ASV),]$adjPvalues < 0.05] <- 1
input.5$`GLM NBZI`[glm_nbzi[order(glm_nbzi$Representative_ASV),]$FDR_BH < 0.05] <- 1
input.5$`limma-voom`[limma_voom[order(limma_voom$Representative_ASV),]$adj.P.Val < 0.05] <- 1
input.5$`SAMseq`[samseq[order(samseq$Representative_ASV),]$q.value... < 5] <- 1
input.5$`LEfSe`[lefse[order(lefse$Representative_ASV),]$FDR_BH != "-"] <- 1
input.5[is.na(input.5)] <- 0
input.5$MRA <- (kw_tss[order(kw_tss$Representative_ASV),]$Case_MRA +
                kw_tss[order(kw_tss$Representative_ASV),]$Control_MRA)/2
input.5$MRAR <- kw_tss[order(kw_tss$Representative_ASV),]$MRAR

##### Replicated signatures unfiltered #####

# Filter for genera common between datasets
input.4.sub <- input.4[rownames(input.4) %in% rownames(input.5),]
input.5.sub  <- input.5[rownames(input.5) %in% rownames(input.4),]

# Add matrices together to see who replicated based on significance and who did not
input.6 <- input.4.sub[order(rownames(input.4.sub)),1:(ncol(input.4.sub)-2)] + input.5.sub[order(rownames(input.5.sub)),1:(ncol(input.5.sub)-2)]

# Differentiate replicated associations with same effect directions vs opposite effect directions
input.6[input.6 == 2 & input.4.sub[order(rownames(input.4.sub)),]$MRAR > 1 & input.5.sub[order(rownames(input.5.sub)),]$MRAR > 1] <- 3
input.6[input.6 == 2 & input.4.sub[order(rownames(input.4.sub)),]$MRAR < 1 & input.5.sub[order(rownames(input.5.sub)),]$MRAR < 1] <- 3

# Convert back to binary matrix
input.6[input.6 != 3] <- 0
input.6[input.6 == 3] <- 1

# Add MRAs and MRARs from both datasets
input.6 <- cbind(input.6, data.frame(MRAR_1=input.4.sub[order(rownames(input.4.sub)),]$MRAR,
                                     MRAR_2=input.5.sub[order(rownames(input.5.sub)),]$MRAR,
                                     MRA_1=input.4.sub[order(rownames(input.4.sub)),]$MRA,
                                     MRA_2=input.5.sub[order(rownames(input.5.sub)),]$MRA))

# Print N replicated
rep2 <- data.frame(Taxa=names(rowSums(input.6[,1:(ncol(input.6)-4)])), 
                   `N methods replicated`=rowSums(input.6[,1:(ncol(input.6)-4)]), check.names=F)
write.table(rep2, "../Script_Output/Joint_Analyses_Output/N_replicated_10%.txt", quote=F, row.names=F, sep="\t")

# Calculate proportion of differential abundance signatures detected by each method out of the total taxa tested
p1 <- round(colSums(input.1[,1:(ncol(input.1)-2)])/nrow(input.1), 2)
p2 <- c(round(colSums(input.2[,1:(ncol(input.2)-2)])/nrow(input.2), 2), SAMseq=NA)
p3 <- c(round(colSums(input.3[,1:(ncol(input.3)-4)])/nrow(input.3), 2), SAMseq=NA)
p4 <- round(colSums(input.4[,1:(ncol(input.4)-2)])/nrow(input.4), 2)
p5 <- round(colSums(input.5[,1:(ncol(input.5)-2)])/nrow(input.5), 2)
p6 <- round(colSums(input.6[,1:(ncol(input.6)-4)])/nrow(input.6), 2)

##### Prepping pairwise concordance data #####

# Read in concordances for all scenarios and set diagonals to NA
d1 <- read.table("../Script_Output/Dataset1_Output/Pairwise_concordances.txt", row.names=1, header=T, sep="\t", check.names=F)
diag(d1[,ncol(d1):1]) <- NA
d1.filt <- read.table("../Script_Output/Dataset1_Output/Pairwise_concordances_10%.txt", row.names=1, header=T, sep="\t", check.names=F)
diag(d1.filt[,ncol(d1.filt):1]) <- NA
d2 <- read.table("../Script_Output/Dataset2_Output/Pairwise_concordances.txt", row.names=1, header=T, sep="\t", check.names=F)
diag(d2[,ncol(d2):1]) <- NA
d2.filt <- read.table("../Script_Output/Dataset2_Output/Pairwise_concordances_10%.txt", row.names=1, header=T, sep="\t", check.names=F)
diag(d2.filt[,ncol(d2.filt):1]) <- NA
d3 <- read.table("../Script_Output/Joint_Analyses_Output/Pairwise_concordances.txt", row.names=1, header=T, sep="\t", check.names=F)
diag(d3[,ncol(d3):1]) <- NA
d3.filt <- read.table("../Script_Output/Joint_Analyses_Output/Pairwise_concordances_10%.txt", row.names=1, header=T, sep="\t", check.names=F)
diag(d3.filt[,ncol(d3.filt):1]) <- NA

# Put NA for SAMseq for dataset 2 and replication aware concordances for unfiltered data
d2 <- cbind(data.frame(SAMseq=rep(NA, nrow(d2))), d2)
d2 <- rbind(d2, rep(NA, ncol(d2)))
rownames(d2)[nrow(d2)] <- "SAMseq"
d3 <- cbind(data.frame(SAMseq=rep(NA, nrow(d3))), d3)
d3 <- rbind(d3, rep(NA, ncol(d3)))
rownames(d3)[nrow(d3)] <- "SAMseq"

# Add descriptive variables
d1$dataset <- "Dataset 1"
d2$dataset <- "Dataset 2"
d3$dataset <- "Replicated"
d1.filt$dataset <- "Dataset 1"
d2.filt$dataset <- "Dataset 2"
d3.filt$dataset <- "Replicated"

# Combine unfiltered and filtered results
unfilt <- rbind(d1, d2, d3)
filt <- rbind(d1.filt, d2.filt, d3.filt)

# Sort data by overall means
means <- sort(apply(unfilt[,1:(ncol(unfilt)-1)], 2, function(x){mean(na.omit(x))}))
unfilt <- unfilt[, match(c(names(means),"dataset"), colnames(unfilt))]

means <- sort(apply(filt[,1:(ncol(filt)-1)], 2, function(x){mean(na.omit(x))}))
filt <- filt[, match(c(names(means),"dataset"), colnames(filt))]

##### Plotting data #####

# Melt data
unfilt <- melt(unfilt, id.vars="dataset")
filt <- melt(filt, id.vars="dataset")

# Add caclulated proportions
unfilt$prop <- NA
for (i in 1:nrow(unfilt[unfilt$dataset == "Dataset 1",])){
  unfilt[unfilt$dataset == "Dataset 1",][i,"prop"] <- p1[as.character(unfilt[unfilt$dataset == "Dataset 1",][i,"variable"])]
}
for (i in 1:nrow(unfilt[unfilt$dataset == "Dataset 2",])){
  unfilt[unfilt$dataset == "Dataset 2",][i,"prop"] <- p2[as.character(unfilt[unfilt$dataset == "Dataset 2",][i,"variable"])]
}
for (i in 1:nrow(unfilt[unfilt$dataset == "Replicated",])){
  unfilt[unfilt$dataset == "Replicated",][i,"prop"] <- p3[as.character(unfilt[unfilt$dataset == "Replicated",][i,"variable"])]
}

filt$prop <- NA
for (i in 1:nrow(filt[filt$dataset == "Dataset 1",])){
  filt[filt$dataset == "Dataset 1",][i,"prop"] <- p4[as.character(filt[filt$dataset == "Dataset 1",][i,"variable"])]
}
for (i in 1:nrow(filt[filt$dataset == "Dataset 2",])){
  filt[filt$dataset == "Dataset 2",][i,"prop"] <- p5[as.character(filt[filt$dataset == "Dataset 2",][i,"variable"])]
}
for (i in 1:nrow(filt[filt$dataset == "Replicated",])){
  filt[filt$dataset == "Replicated",][i,"prop"] <- p6[as.character(filt[filt$dataset == "Replicated",][i,"variable"])]
}

# Double data making one half data for boxplots, and the other half for bar plot
unfilt <- rbind(data.frame(unfilt, plot="1"), data.frame(unfilt, plot="2"))
unfilt$value[unfilt$plot == 2] <- NA

filt <- rbind(data.frame(filt, plot="1"), data.frame(filt, plot="2"))
filt$value[filt$plot == 2] <- NA

# Add mean concordance to plot for boxplots of concordances
unfilt$mean <- NA
unfilt$mean[unfilt$dataset == "Dataset 1" & unfilt$plot == 1] <- mean(na.omit(unfilt$value[unfilt$dataset == "Dataset 1"]))
unfilt$mean[unfilt$dataset == "Dataset 2" & unfilt$plot == 1] <- mean(na.omit(unfilt$value[unfilt$dataset == "Dataset 2"]))
unfilt$mean[unfilt$dataset == "Replicated" & unfilt$plot == 1] <- mean(na.omit(unfilt$value[unfilt$dataset == "Replicated"]))

filt$mean <- NA
filt$mean[filt$dataset == "Dataset 1" & filt$plot == 1] <- mean(na.omit(filt$value[filt$dataset == "Dataset 1"]))
filt$mean[filt$dataset == "Dataset 2" & filt$plot == 1] <- mean(na.omit(filt$value[filt$dataset == "Dataset 2"]))
filt$mean[filt$dataset == "Replicated" & filt$plot == 1] <- mean(na.omit(filt$value[filt$dataset == "Replicated"]))

# Add color for plotting
col.pal <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#332288","#117733",
             "#964B00","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255","#648FFF","#785EF0","#DC267F","#FE6100",
             "#FFB000","#E1BE6A","#999999")
set.seed(1234)
names(col.pal) <- sample(levels(unfilt$variable))

unfilt$color <- NA
for (i in 1:nrow(unfilt)){
  unfilt[i,"color"] <- col.pal[as.character(unfilt[i,"variable"])]
}

filt$color <- NA
for (i in 1:nrow(filt)){
  filt[i,"color"] <- col.pal[as.character(filt[i,"variable"])]
}

# Get correlation labels for dot plots
d1.cor <- cor.test(unfilt$prop[unfilt$dataset == "Dataset 1"], unfilt$value[unfilt$dataset == "Dataset 1"])
d2.cor <- cor.test(unfilt$prop[unfilt$dataset == "Dataset 2"], unfilt$value[unfilt$dataset == "Dataset 2"])
d3.cor <- cor.test(unfilt$prop[unfilt$dataset == "Replicated"], unfilt$value[unfilt$dataset == "Replicated"])
unfilt.cor.lab <- data.frame(label = c(paste("Pearson's r [95% CI] = ",round(d1.cor$estimate,2)," [",round(d1.cor$conf.int[1],2),"; ",round(d1.cor$conf.int[2],2),"]",", ",
                                             "P = ",signif(d1.cor$p.value,0),sep=""),
                                       paste("Pearson's r [95% CI] = ",round(d2.cor$estimate,2)," [",round(d2.cor$conf.int[1],2),"; ",round(d2.cor$conf.int[2],2),"]",", ",
                                             "P = ",signif(d2.cor$p.value,0),sep=""),
                                       paste("Pearson's r [95% CI] = ",round(d3.cor$estimate,2)," [",round(d3.cor$conf.int[1],2),"; ",round(d3.cor$conf.int[2],2),"]",", ",
                                             "P = ",signif(d3.cor$p.value,0),sep="")),
                             dataset = c("Dataset 1", "Dataset 2", "Replicated"))

d1.cor <- cor.test(filt$prop[filt$dataset == "Dataset 1"], filt$value[filt$dataset == "Dataset 1"])
d2.cor <- cor.test(filt$prop[filt$dataset == "Dataset 2"], filt$value[filt$dataset == "Dataset 2"])
d3.cor <- cor.test(filt$prop[filt$dataset == "Replicated"], filt$value[filt$dataset == "Replicated"])
filt.cor.lab <- data.frame(label = c(paste("Pearson's r [95% CI] = ",round(d1.cor$estimate,2)," [",round(d1.cor$conf.int[1],2),"; ",round(d1.cor$conf.int[2],2),"]",", ",
                                             "P = ",signif(d1.cor$p.value,0),sep=""),
                                       paste("Pearson's r [95% CI] = ",round(d2.cor$estimate,2)," [",round(d2.cor$conf.int[1],2),"; ",round(d2.cor$conf.int[2],2),"]",", ",
                                             "P = ",signif(d2.cor$p.value,0),sep=""),
                                       paste("Pearson's r [95% CI] = ",round(d3.cor$estimate,2)," [",round(d3.cor$conf.int[1],2),"; ",round(d3.cor$conf.int[2],2),"]",", ",
                                             "P = ",signif(d3.cor$p.value,0),sep="")),
                             dataset = c("Dataset 1", "Dataset 2", "Replicated"))

# For filtered data, add difference in mean concordances between unfiltered and filtered
mean1 <- aggregate(value ~ variable, mean, data=filt[filt$dataset == "Dataset 1",])
mean2 <- aggregate(value ~ variable, mean, data=unfilt[unfilt$dataset == "Dataset 1",])
d1.diff <- round(mean1$value - mean2[match(mean1$variable, mean2$variable),]$value, 2)
d1.diff <- data.frame(dataset="Dataset 1", plot=1, variable=mean1$variable, diff=d1.diff)
d1.diff$color <- ifelse(d1.diff$diff < 0, "2", "1")
d1.diff$diff <- ifelse(d1.diff$diff > 0, paste("+",d1.diff$diff,sep=""), d1.diff$diff)

mean1 <- aggregate(value ~ variable, mean, data=filt[filt$dataset == "Dataset 2",])
mean2 <- aggregate(value ~ variable, mean, data=unfilt[unfilt$dataset == "Dataset 2",])
d2.diff <- round(mean1$value - mean2[match(mean1$variable, mean2$variable),]$value, 2)
d2.diff <- data.frame(dataset="Dataset 2", plot=1, variable=mean1$variable, diff=d2.diff)
d2.diff$color <- ifelse(d2.diff$diff < 0, "2", "1")
d2.diff$diff <- ifelse(d2.diff$diff > 0, paste("+",d2.diff$diff,sep=""), d2.diff$diff)

mean1 <- aggregate(value ~ variable, mean, data=filt[filt$dataset == "Replicated",])
mean2 <- aggregate(value ~ variable, mean, data=unfilt[unfilt$dataset == "Replicated",])
d3.diff <- round(mean1$value - mean2[match(mean1$variable, mean2$variable),]$value, 2)
d3.diff <- data.frame(dataset="Replicated", plot=1, variable=mean1$variable, diff=d3.diff)
d3.diff$color <- ifelse(d3.diff$diff < 0, "2", "1")
d3.diff$diff <- ifelse(d3.diff$diff > 0, paste("+",d3.diff$diff,sep=""), d3.diff$diff)

filt.conc.diff.lab <- rbind(d1.diff, d2.diff, d3.diff)

# For filtered data, add difference in proportions between unfiltered and filtered
mean1 <- aggregate(prop ~ variable, mean, data=filt[filt$dataset == "Dataset 1",])
mean2 <- aggregate(prop ~ variable, mean, data=unfilt[unfilt$dataset == "Dataset 1",])
d1.diff <- round(mean1$prop - mean2[match(mean1$variable, mean2$variable),]$prop, 2)
d1.diff <- data.frame(dataset="Dataset 1", plot=2, variable=mean1$variable, diff=d1.diff)
d1.diff$color <- ifelse(d1.diff$diff < 0, "2", "1")
d1.diff$diff <- ifelse(d1.diff$diff > 0, paste("+",d1.diff$diff,sep=""), d1.diff$diff)

mean1 <- aggregate(prop ~ variable, mean, data=filt[filt$dataset == "Dataset 2",])
mean2 <- aggregate(prop ~ variable, mean, data=unfilt[unfilt$dataset == "Dataset 2",])
d2.diff <- round(mean1$prop - mean2[match(mean1$variable, mean2$variable),]$prop, 2)
d2.diff <- data.frame(dataset="Dataset 2", plot=2, variable=mean1$variable, diff=d2.diff)
d2.diff$color <- ifelse(d2.diff$diff < 0, "2", "1")
d2.diff$diff <- ifelse(d2.diff$diff > 0, paste("+",d2.diff$diff,sep=""), d2.diff$diff)

mean1 <- aggregate(prop ~ variable, mean, data=filt[filt$dataset == "Replicated",])
mean2 <- aggregate(prop ~ variable, mean, data=unfilt[unfilt$dataset == "Replicated",])
d3.diff <- round(mean1$prop - mean2[match(mean1$variable, mean2$variable),]$prop, 2)
d3.diff <- data.frame(dataset="Replicated", plot=2, variable=mean1$variable, diff=d3.diff)
d3.diff$color <- ifelse(d3.diff$diff < 0, "2", "1")
d3.diff$diff <- ifelse(d3.diff$diff > 0, paste("+",d3.diff$diff,sep=""), d3.diff$diff)

filt.prop.diff.lab <- rbind(d1.diff, d2.diff, d3.diff)

# Plot data for results with unfiltered data
pdf("../Script_Output/Joint_Analyses_Output/Concordance_and_proportion_DA_summary.pdf", onefile=F, height=10, width=8)

g1 <- ggplot(data = unfilt[unfilt$plot==1,], aes(x = variable, y = value)) + 
        geom_boxplot(fill="grey", notch=F, outlier.size=-1) + 
        geom_dotplot(inherit.aes=F, data = unfilt[unfilt$plot==1,], aes(x = variable, y = value), binaxis = "y", 
                     stackdir = "center", binwidth = 0.03, position = position_dodge(0.75), dotsize = 0.5) +
        stat_summary(fun=mean, geom="point", shape=1, size=1, color="red") +
        geom_hline(data=unfilt[unfilt$plot==1,], aes(yintercept=mean), color="red", linetype="solid", size=0.25) +
        geom_bar(inherit.aes=F, data=unique(unfilt[unfilt$plot==2,]), aes(x = variable, y = prop, fill=variable), stat="identity") +
        facet_nested(dataset + plot ~ ., labeller=labeller(plot=c("1"="Concordances", "2"="Proportion DA"))) +
        scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1)) +
        scale_x_discrete(position="top") +
        scale_fill_manual(values=unique(unfilt[unfilt$plot==1,]$color)) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 25, hjust = 0, size = 5), axis.text.x.top = element_text(vjust = 0.05)) + 
        theme(axis.text.y = element_text(size = 6)) + 
        theme(strip.text = element_text(margin = margin(0,0.01,0,0.01, "in"), size = 7)) +
        theme(plot.margin = unit(c(-0.75,1,1.25,0), "lines")) +
        guides(fill=F) + 
        labs(x="", y="") 

g2 <- ggplot(data = unfilt[unfilt$plot==1,], aes(x = prop, y = value)) + 
        stat_smooth(method="lm", n=nrow(unfilt[unfilt$plot==1,]), color="black", alpha=0.8) +
        geom_point(aes(color = variable), size=1) +
        facet_nested(dataset ~ .) +
        scale_color_manual(values=unique(unfilt[unfilt$plot==1,]$color)) +
        scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
        coord_cartesian(xlim=c(0,0.8), ylim=c(0,1)) +
        geom_text(data=unfilt.cor.lab, aes(x=0.5, y=1, label=label), size=2) +
        scale_x_continuous(position="top") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 6), axis.text.x = element_text(size = 6)) + 
        theme(axis.title.y = element_text(size = 6), axis.text.y = element_text(size = 6)) + 
        theme(strip.text = element_text(margin = margin(0,0.01,0,0.01, "in"), size = 7)) +
        theme(plot.margin = unit(c(0.7,0.5,1.25,-0.5), "lines")) +
        theme(legend.position="right", legend.title=element_blank(), legend.text = element_text(size=6), legend.spacing.x = unit(0.01, "in"), 
              legend.box.margin=margin(-10,-10,-10,-10)) +
        guides(color=guide_legend(reverse=T, ncol=1, keyheight = 1, keywidth = 0.5, override.aes = list(size=2))) +
        labs(x="Proportion DA", y="Concordances") 

ggarrange(g1, g2, ncol=2, labels=c("A","B"))

dev.off()

# Plot data for results with filtered data
pdf(onefile=F,height=10, width=8)

g1 <- ggplot(data = filt[filt$plot==1,], aes(x = variable, y = value)) + 
        geom_boxplot(fill="grey", notch=F, outlier.size=-1) + 
        geom_dotplot(inherit.aes=F, data = filt[filt$plot==1,], aes(x = variable, y = value), binaxis = "y", 
                     stackdir = "center", binwidth = 0.03, position = position_dodge(0.75), dotsize = 0.5) +
        stat_summary(fun=mean, geom="point", shape=1, size=1, color="red") +
        geom_hline(data=filt[filt$plot==1,], aes(yintercept=mean), color="red", linetype="solid", size=0.25) +
        geom_bar(inherit.aes=F, data=unique(filt[filt$plot==2,]), aes(x = variable, y = prop, fill=variable), stat="identity") +
        geom_text(inherit.aes=F, data=filt.conc.diff.lab, aes(x=variable, y=1.1, label=diff, color=color), angle=55, size=2) +
        geom_text(inherit.aes=F, data=filt.prop.diff.lab, aes(x=variable, y=1.1, label=diff, color=color), angle=55, size=2) +
        facet_nested(dataset + plot ~ ., labeller=labeller(plot=c("1"="Concordances", "2"="Proportion DA"))) +
        scale_y_continuous(limits=c(0,1.15), breaks=c(0,0.25,0.5,0.75,1)) +
        scale_x_discrete(position="top") +
        scale_fill_manual(values=unique(filt[filt$plot==1,]$color)) +
        scale_color_manual(values=c("darkgreen","darkred")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 25, hjust = 0, size = 4.5), axis.text.x.top = element_text(vjust = 0.05)) + 
        theme(axis.text.y = element_text(size = 5.5)) + 
        theme(strip.text = element_text(margin = margin(0,0.01,0,0.01, "in"), size = 6.5)) +
        theme(plot.margin = unit(c(-0.75,1,1.25,0), "lines")) +
        guides(color=F, fill=F) + 
        labs(x="", y="") 

g2 <- ggplot(data = filt[filt$plot==1,], aes(x = prop, y = value)) + 
        stat_smooth(method="lm", n=nrow(filt[filt$plot==1,]), color="black", alpha=0.8) +
        geom_point(aes(color = variable), size=1) +
        facet_nested(dataset ~ .) +
        scale_color_manual(values=unique(filt[filt$plot==1,]$color)) +
        scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) +
        coord_cartesian(xlim=c(0,0.8), ylim=c(0,1.1)) +
        geom_text(data=filt.cor.lab, aes(x=0.5, y=1.1, label=label), size=2) +
        scale_x_continuous(position="top") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5.5), axis.text.x = element_text(size = 5.5)) + 
        theme(axis.title.y = element_text(size = 5.5), axis.text.y = element_text(size = 5.5)) + 
        theme(strip.text = element_text(margin = margin(0,0.01,0,0.01, "in"), size = 6.5)) +
        theme(plot.margin = unit(c(0.7,0.5,1.25,-0.5), "lines")) +
        theme(legend.position="right", legend.title=element_blank(), legend.text = element_text(size=5.5), legend.spacing.x = unit(0.01, "in"), 
              legend.box.margin=margin(-10,-10,-10,-10)) +
        guides(color=guide_legend(reverse=T, ncol=1, keyheight = 1, keywidth = 0.5, override.aes = list(size=2))) +
        labs(x="Proportion DA", y="Concordances") 

ggarrange(g1, g2, ncol=2, labels=c("A","B"))

dev.off()
system("mv Rplot*.pdf ../Script_Output/Joint_Analyses_Output/Concordance_and_proportion_DA_summary_10%.pdf")

###########################################################################
##### Create plot summarizing the distribution of taxa detected as    #####
##### as a function of mean relative abundance and fold change in PD  #####
###########################################################################

# Prep binary taxa by method matrix data for plotting
input.1.plot <- melt(data.frame(input.1[,1:(ncol(input.1)-2)], dataset="Dataset 1", check.names=F), id.vars="dataset")
input.1.plot$MRA <- rep(input.1$MRA, 23)
input.1.plot$MRAR <- rep(input.1$MRAR, 23)
input.1.plot <- input.1.plot[is.finite(input.1.plot$MRAR),]
input.1.plot <- input.1.plot[input.1.plot$MRAR > 0,]

input.2.plot <- melt(data.frame(input.2[,1:(ncol(input.2)-2)], SAMseq=NA, dataset="Dataset 2", check.names=F), id.vars="dataset")
input.2.plot$MRA <- rep(input.2$MRA, 23)
input.2.plot$MRAR <- rep(input.2$MRAR, 23)
input.2.plot <- input.2.plot[is.finite(input.2.plot$MRAR),]
input.2.plot <- input.2.plot[input.2.plot$MRAR > 0,]

input.3.plot <- melt(data.frame(input.3[,1:(ncol(input.3)-4)], SAMseq=NA, dataset="Replicated", check.names=F), id.vars="dataset")
input.3.plot$MRA <- rep((input.3$MRA_1+input.3$MRA_2)/2, 23)
input.3.plot$MRAR <- rep((input.3$MRAR_1+input.3$MRAR_2)/2, 23)
input.3.plot <- input.3.plot[is.finite(input.3.plot$MRAR),]
input.3.plot <- input.3.plot[input.3.plot$MRAR > 0,]

input.4.plot <- melt(data.frame(input.4[,1:(ncol(input.4)-2)], dataset="Dataset 1", check.names=F), id.vars="dataset")
input.4.plot$MRA <- rep(input.4$MRA, 23)
input.4.plot$MRAR <- rep(input.4$MRAR, 23)
input.4.plot <- input.4.plot[is.finite(input.4.plot$MRAR),]
input.4.plot <- input.4.plot[input.4.plot$MRAR > 0,]

input.5.plot <- melt(data.frame(input.5[,1:(ncol(input.5)-2)], dataset="Dataset 2", check.names=F), id.vars="dataset")
input.5.plot$MRA <- rep(input.5$MRA, 23)
input.5.plot$MRAR <- rep(input.5$MRAR, 23)
input.5.plot <- input.5.plot[is.finite(input.5.plot$MRAR),]
input.5.plot <- input.5.plot[input.5.plot$MRAR > 0,]

input.6.plot <- melt(data.frame(input.6[,1:(ncol(input.6)-4)], dataset="Replicated", check.names=F), id.vars="dataset")
input.6.plot$MRA <- rep((input.6$MRA_1+input.6$MRA_2)/2, 23)
input.6.plot$MRAR <- rep((input.6$MRAR_1+input.6$MRAR_2)/2, 23)
input.6.plot <- input.6.plot[is.finite(input.6.plot$MRAR),]
input.6.plot <- input.6.plot[input.6.plot$MRAR > 0,]

# Combine results from unfiltered and filtered
unfilt2 <- rbind(input.1.plot, input.2.plot, input.3.plot)
filt2 <- rbind(input.4.plot, input.5.plot, input.6.plot)

# Remove NA entries
unfilt2 <- unfilt2[!is.na(unfilt2$value),]
filt2 <- filt2[!is.na(filt2$value),]

# Order based on previous figure
unfilt2$variable <- factor(unfilt2$variable, levels=rev(levels(unfilt$variable)))
filt2$variable <- factor(filt2$variable, levels=rev(levels(filt$variable)))

# Add median mean relative abundance for plotting
unfilt2$median <- NA
unfilt2$median[unfilt2$dataset == "Dataset 1"] <- median(na.omit(unfilt2$MRA[unfilt2$dataset == "Dataset 1"]))
unfilt2$median[unfilt2$dataset == "Dataset 2"] <- median(na.omit(unfilt2$MRA[unfilt2$dataset == "Dataset 2"]))
unfilt2$median[unfilt2$dataset == "Replicated"] <- median(na.omit(unfilt2$MRA[unfilt2$dataset == "Replicated"]))

filt2$median <- NA
filt2$median[filt2$dataset == "Dataset 1"] <- median(na.omit(filt2$MRA[filt2$dataset == "Dataset 1"]))
filt2$median[filt2$dataset == "Dataset 2"] <- median(na.omit(filt2$MRA[filt2$dataset == "Dataset 2"]))
filt2$median[filt2$dataset == "Replicated"] <- median(na.omit(filt2$MRA[filt2$dataset == "Replicated"]))

# Add variable that tags what MRA is above or below median
unfilt2$prev <- ifelse(unfilt2$MRA > unfilt2$median, 1, 0)
filt2$prev <- ifelse(filt2$MRA > filt2$median, 1, 0)

# Add variable that tags what MRAR is within or outside a log2 fold change of 0.4 (~1.3x change)
unfilt2$fc_area <- ifelse(log2(unfilt2$MRAR) < -0.4 | log2(unfilt2$MRAR) > 0.4, 1, 0)
filt2$fc_area <- ifelse(log2(filt2$MRAR) < -0.4 | log2(filt2$MRAR) > 0.4, 1, 0)

# Perform Fisher exact test for enrichment of certain taxa as DA signatures based on MRA and log2(MRAR)
unfilt_lab <- data.frame()
for (i in 1:length(levels(unfilt2$variable))){
  curr.method <- levels(unfilt2$variable)[i]
  for (j in 1:length(levels(as.factor(unfilt2$dataset)))){
    curr.dataset <- levels(as.factor(unfilt2$dataset))[j]
    if ((curr.method == "SAMseq" & curr.dataset == "Dataset 2") | (curr.method == "SAMseq" & curr.dataset == "Replicated")){
      next
    }else{
      sub.data <- unfilt2[unfilt2$variable == curr.method & unfilt2$dataset == curr.dataset,]
      f.out.prev <- fisher.test(table(sub.data$prev, sub.data$value))
      f.out.fc <- fisher.test(table(sub.data$fc_area, sub.data$value))
      unfilt_lab <- rbind(unfilt_lab, data.frame(dataset=curr.dataset, variable=curr.method,
                                                 prev_lab=paste("MRA: OR=",round(f.out.prev$estimate, 2),", P=",signif(f.out.prev$p.value, 0),sep=""),
                                                 fc_lab=paste("FC: OR=",round(f.out.fc$estimate, 2),", P=",signif(f.out.fc$p.value, 0),sep="")))
    }
  }
}
unfilt_lab$variable <- factor(unfilt_lab$variable, levels=levels(unfilt2$variable))

filt_lab <- data.frame()
for (i in 1:length(levels(filt2$variable))){
  curr.method <- levels(filt2$variable)[i]
  for (j in 1:length(levels(as.factor(filt2$dataset)))){
    curr.dataset <- levels(as.factor(filt2$dataset))[j]
    sub.data <- filt2[filt2$variable == curr.method & filt2$dataset == curr.dataset,]
    f.out.prev <- fisher.test(table(sub.data$prev, sub.data$value))
    f.out.fc <- fisher.test(table(sub.data$fc_area, sub.data$value))
    filt_lab <- rbind(filt_lab, data.frame(dataset=curr.dataset, variable=curr.method,
                                           prev_lab=paste("MRA: OR=",round(f.out.prev$estimate, 2),", P=",signif(f.out.prev$p.value, 0),sep=""),
                                           fc_lab=paste("FC: OR=",round(f.out.fc$estimate, 2),", P=",signif(f.out.fc$p.value, 0),sep="")))
  }
}
filt_lab$variable <- factor(filt_lab$variable, levels=levels(filt2$variable))

# Plot data for results with unfiltered data
pdf("../Script_Output/Joint_Analyses_Output/MRAR_vs_MRA.pdf", onefile=F, height=11, width=5)

g1 <- ggplot(data = unfilt2[order(unfilt2$value),], aes(x = log2(MRAR), y = MRA, color=as.character(value))) + 
        geom_point(size=0.25) +
        geom_vline(data=unfilt2[order(unfilt2$value),], aes(xintercept=0.4), color="red", linetype="dashed", size=0.25) +
        geom_vline(data=unfilt2[order(unfilt2$value),], aes(xintercept=-0.4), color="red", linetype="dashed", size=0.25) +
        geom_hline(data=unfilt2[order(unfilt2$value),], aes(yintercept=median), color="red", linetype="dashed", size=0.25) +
        geom_text(inherit.aes=F, data=unfilt_lab, aes(x=9.5, y=0.1, label=prev_lab), size=1.5) +
        geom_text(inherit.aes=F, data=unfilt_lab, aes(x=10, y=2e-7, label=fc_lab), size=1.5) +
        facet_grid(variable ~ dataset) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
        scale_color_manual(values=c("#999999", "#000000")) +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5), axis.text.x = element_text(size = 4)) + 
        theme(axis.title.y = element_text(size = 5), axis.text.y = element_text(size = 4)) + 
        theme(strip.text = element_text(margin = margin(0,0.01,0,0.01, "in"), size = 3.5)) +
        theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
        theme(panel.spacing = unit(0.1, "lines")) +
        guides(color=F) + 
        labs(x="log2 fold change in PD", y="Mean relative abundances (log scale)") 

g1 <- ggplot_gtable(ggplot_build(g1))
stripr <- which(grepl('strip-r', g1$layout$name))
fills <- c(rep("#CCF5F7",10), rep("grey80", 4), rep("#F7D2CC", 9))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g1)

dev.off()

# Plot data for results with filtered data
pdf(onefile=F, height=11, width=5)

g1 <- ggplot(data = filt2[order(filt2$value),], aes(x = log2(MRAR), y = MRA, color=as.character(value))) + 
        geom_point(size=0.25) +
        geom_vline(data=filt2[order(filt2$value),], aes(xintercept=0.4), color="red", linetype="dashed", size=0.25) +
        geom_vline(data=filt2[order(filt2$value),], aes(xintercept=-0.4), color="red", linetype="dashed", size=0.25) +
        geom_hline(data=filt2[order(filt2$value),], aes(yintercept=median), color="red", linetype="dashed", size=0.25) +
        geom_text(inherit.aes=F, data=filt_lab, aes(x=6, y=0.1, label=prev_lab), size=1.5) +
        geom_text(inherit.aes=F, data=filt_lab, aes(x=6, y=0.025, label=fc_lab), size=1.5) +
        facet_grid(variable ~ dataset) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
        scale_color_manual(values=c("#999999", "#000000")) +
        theme_bw() +
        theme(axis.title.x = element_text(size = 5), axis.text.x = element_text(size = 4)) + 
        theme(axis.title.y = element_text(size = 5), axis.text.y = element_text(size = 4)) + 
        theme(strip.text = element_text(margin = margin(0,0.01,0,0.01, "in"), size = 3.5)) +
        theme(plot.margin = unit(c(0.5,0.5,0,0), "lines")) +
        theme(panel.spacing = unit(0.1, "lines")) +
        guides(color=F) + 
        labs(x="log2 fold change in PD", y="Mean relative abundances (log scale)") 

g1 <- ggplot_gtable(ggplot_build(g1))
stripr <- which(grepl('strip-r', g1$layout$name))
fills <- c(rep("#CCF5F7",11), "grey80", "grey80", "grey80", rep("#F7D2CC", 9))
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
  g1$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g1)

dev.off()
system("mv Rplot*.pdf ../Script_Output/Joint_Analyses_Output/MRAR_vs_MRA_10%.pdf")

######################################################################
##### Create taxa by method heatmap with hierarchical clustering #####
##### using results from performing DA testing on filtered data  #####
######################################################################

input <- input.6
                                            
# Create colors for control relative abundances
relabun.colors <- colorRampPalette(c("white", "#00810A"))
input <- input[order(input$MRA_1),]
input$MRA_1_col <- relabun.colors(dim(input)[1])
input <- input[order(input$MRA_2),]
input$MRA_2_col <- relabun.colors(dim(input)[1])

# Create colors for fold changes
down.colors <- colorRampPalette(c("#BF3019","white"))
up.colors <- colorRampPalette(c("white","#0066BA"))
input <- input[order(input$MRAR_1),]
input$MRAR_1_col[input$MRAR_1 < 1] <- down.colors(dim(input[input$MRAR_1 < 1,])[1])
input$MRAR_1_col[input$MRAR_1 > 1] <- up.colors(dim(input[input$MRAR_1 > 1,])[1])
input <- input[order(input$MRAR_2),]
input$MRAR_2_col[input$MRAR_2 < 1] <- down.colors(dim(input[input$MRAR_2 < 1,])[1])
input$MRAR_2_col[input$MRAR_2 > 1] <- up.colors(dim(input[input$MRAR_2 > 1,])[1])

# Filter out taxa who were not replicated by at least 4 methods
input <- input[rowSums(input[,1:(dim(input)[2]-8)]) >= 1,]

# Create more readable taxa labels
taxonomy <- strsplit(rownames(input), ";")
tax.lab <- c()
for (i in 1:length(taxonomy)){
  curr.taxa <- taxonomy[[i]]
  if (sum(curr.taxa == "NA") == 0){
    tax.lab <- append(tax.lab, curr.taxa[6])
  }else if (sum(curr.taxa == "NA") == 1){
    tax.lab <- append(tax.lab, paste(curr.taxa[5],"_unclass",sep=""))
  }else if (sum(curr.taxa == "NA") == 2){
    tax.lab <- append(tax.lab, paste(curr.taxa[4],"_unclass",sep=""))
  }else if (sum(curr.taxa == "NA") == 3){
    tax.lab <- append(tax.lab, paste(curr.taxa[3],"_unclass",sep=""))
  }else if (sum(curr.taxa == "NA") == 4){
    tax.lab <- append(tax.lab, paste(curr.taxa[2],"_unclass",sep=""))
  }else if (sum(curr.taxa == "NA") == 5){
    tax.lab <- append(tax.lab, paste(curr.taxa[1],"_unclass",sep=""))
  }else{
    tax.lab <- append(tax.lab, paste("unclassified",sep=""))
  }
}
rownames(input) <- tax.lab

# Plot heatmap using euclidean distance and diana clustering
pdf("../Script_Output/Joint_Analyses_Output/GenusxMethod_heatmap.pdf", height=6, width=6)

hm <- heatmap.3(x=as.matrix(input[,1:(dim(input)[2]-8)]), distfun=function(x)dist(x, method="binary"), hclustfun=function(x)as.hclust(diana(x)), 
                dendrogram="both", col=c("grey90","black"),
                RowSideColors=t(matrix(c(input$MRAR_1_col, input$MRAR_2_col, input$MRA_1_col, input$MRA_2_col), 
                                       nrow=dim(input)[1], ncol=4)), 
                cexRow=0.5, cexCol=0.75, RowSideColorsSize=3, key=F, margins=c(8,8))

lgd = rep(NA, 11)
lgd[c(1,6,11)] = c(7.6,1,0.4)
legend(x=0, y=1, legend = lgd, fill = colorRampPalette(colors = c('#0066BA','grey90','#BF3019'))(11),
       border = F, y.intersp = 0.5, cex = 0.5)
                                           
lgd = rep(NA, 11)
lgd[c(1,11)] = c(0.22,1e-4)
legend(x=0.1, y=1, legend = lgd, fill = colorRampPalette(colors = c('#00810A','grey90'))(11),
       border = F, y.intersp = 0.5, cex = 0.5)
                                            
legend("topright", legend=c("Replicated","Not replicated"), fill=c("black","grey90"), border=T, bty="n", y.intersp = 0.8, cex=0.5)
                                            
mtext("MRAR_1", cex=0.4, side=2, line=-4, at=-0.04) #If anything about the input data or heatmap parameters change these will probably have to be manually changed by trial and error

mtext("MRAR_2", cex=0.4, side=2, line=-4.5, at=-0.04) #If anything about the input data or heatmap parameters change these will probably have to be manually changed by trial and error

mtext("MRA_1", cex=0.4, side=2, line=-5, at=-0.035) #If anything about the input data or heatmap parameters change these will probably have to be manually changed by trial and error

mtext("MRA_2", cex=0.4, side=2, line=-5.5, at=-0.035) #If anything about the input data or heatmap parameters change these will probably have to be manually changed by trial and error

dev.off()

# Label what genera were replicated or not for table
input.tab <- input[,1:(dim(input)[2]-8)]
input.tab[input.tab == 0] <- "-"
input.tab[input.tab == 1] <- "Replicated"
input.tab <- cbind(data.frame(MRAR_1=input$MRAR_1, MRAR_2=input$MRAR_2, MRA_1=input$MRA_1, MRA_2=input$MRA_2), input.tab)

# Order columns by leafs of column dendrogram to match heatmap order
input.tab <- input.tab[,c(1:4, match(labels(hm$colDendrogram), colnames(input.tab)))]

# Order rows by leafs of row dendrogram to match heatmap order
input.tab <- input.tab[match(rev(labels(hm$rowDendrogram)), rownames(input.tab)),]

# Write out table describing heatmap values and what genera were detected by what methods in what dataset(s)
write.table(rownames_to_column(input.tab, "Taxa"), "../Script_Output/Joint_Analyses_Output/GenusxMethod_heatmap_table.txt", row.names=F, quote=F, sep="\t")
