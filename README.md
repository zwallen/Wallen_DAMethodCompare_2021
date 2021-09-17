This github repository houses source code used in the manuscript
**Wallen, Z.D. Comparison study of differential abundance testing methods using two large Parkinson disease gut microbiome datasets derived from 16S amplicon sequencing. BMC Bioinformatics 22, 265 (2021). https://doi.org/10.1186/s12859-021-04193-6.**

The following README gives an overview of the overall structure of the repository, and important notes on how to run the scripts.

## Directory tree for repository
```
Wallen_DAMethodCompare_2021
|
|-- BMC_Bioinformatics_Published_Material -- contains manuscript (pdf) and supplementary material (pdf, xlsx) published in BMC Bioinformatics
|
|-- DAreport -- contains an R function (DAreport) and program (run_DAreport.R) that can be used to perform the type of comparisons detailed
|               in the manuscript on a user supplied abundance table. see README for further details.
|
|-- Dataset1_Scripts -- houses scripts used to run differential abundance tests and calculate pairwise concordances for Dataset 1
|   |
|   |-- ALDEx2 -- contains scripts used to run ALDEx2
|   |
|   |-- ANCOM-BC -- contains scripts used to run ANCOM-BC
|   |
|   |-- ANCOM -- contains scripts used to run ANCOM
|   |
|   |-- DESeq2 -- contains scripts used to run DESeq2
|   |
|   |-- GLM -- contains scripts used to run GLM w/ TSS, CLR, and rCLR transform
|   |
|   |-- GLM_NBZI -- contains scripts used to run zero-inflated negative binomial GLM
|   |
|   |-- KW -- contains scripts used to run Kruskal-Wallis rank-sum test w/ TSS, CLR, and rCLR transform
|   |
|   |-- LEfSe -- contains scripts used to run LEfSe
|   |
|   |-- SAMseq -- contains scripts used to run SAMseq
|   |
|   |-- baySeq -- contains scripts used to run baySeq
|   |
|   |-- edgeR_RLE -- contains scripts used to run edgeR w/ RLE
|   |
|   |-- edgeR_TMM -- contains scripts used to run edgeR w/ TMM
|   |
|   |-- fitFeatureModel -- contains scripts used to run fitFeatureModel
|   |
|   |-- fiZIG -- contains scripts used to run fitZIG
|   |
|   |-- limma_voom -- contains scripts used to run limma voom
|   |
|   |-- t_test -- contains scripts used to run Welch's t-test w/ TSS, CLR, and rCLR transform
|
|-- Dataset2_Scripts -- houses scripts used to run differential abundance tests and calculate pairwise concordances for Dataset 2
|   |
|   |-- ALDEx2 -- contains scripts used to run ALDEx2
|   |
|   |-- ANCOM-BC -- contains scripts used to run ANCOM-BC
|   |
|   |-- ANCOM -- contains scripts used to run ANCOM
|   |
|   |-- DESeq2 -- contains scripts used to run DESeq2
|   |
|   |-- GLM -- contains scripts used to run GLM w/ TSS, CLR, and rCLR transform
|   |
|   |-- GLM_NBZI -- contains scripts used to run zero-inflated negative binomial GLM
|   |
|   |-- KW -- contains scripts used to run Kruskal-Wallis rank-sum test w/ TSS, CLR, and rCLR transform
|   |
|   |-- LEfSe -- contains scripts used to run LEfSe
|   |
|   |-- SAMseq -- contains scripts used to run SAMseq
|   |
|   |-- baySeq -- contains scripts used to run baySeq
|   |
|   |-- edgeR_RLE -- contains scripts used to run edgeR w/ RLE
|   |
|   |-- edgeR_TMM -- contains scripts used to run edgeR w/ TMM
|   |
|   |-- fitFeatureModel -- contains scripts used to run fitFeatureModel
|   |
|   |-- fiZIG -- contains scripts used to run fitZIG
|   |
|   |-- limma_voom -- contains scripts used to run limma voom
|   |
|   |-- t_test -- contains scripts used to run Welch's t-test w/ TSS, CLR, and rCLR transform
|
|-- Joint_Analyses_Scripts -- houses scripts used to calculate concordances for replicated DA signatures and create figures
|  
|-- PhyloseqObjects -- phyloseq objects used as input to differential abundance testing
|   |
|   |-- Dataset1 -- phyloseq object for Dataset 1
|   |
|   |-- Dataset2 -- phyloseq object for Dataset 2
|
|-- Script_Output -- directory to house output from scripts in Dataset1_Scripts, Dataset2_Scripts, and Joint_Analyses_Scripts directories
|   |
|   |-- Dataset1_Output -- output from scripts in Dataset1_Scripts directory
|   |   |
|   |   |-- ErrorOut -- stderr output from running scripts
|   |   |
|   |   |-- Output -- stdout output from running scripts
|   |
|   |-- Dataset2_Output -- output from scripts in Dataset2_Scripts directory
|   |   |
|   |   |-- ErrorOut -- stderr output from running scripts
|   |   |
|   |   |-- Output -- stdout output from running scripts
|   |
|   |-- Joint_Analyses_Output -- output from scripts in Joint_Analyses_Scripts directory
|       |
|       |-- ErrorOut -- stderr output from running scripts
|       |
|       |-- Output -- stdout output from running scripts
|
|-- Support_Files -- source code called upon by certain scripts including updated ANCOM code, heatmap.3 code, and function for rCLR transform
```

## Important notes about this repository

#### The repository is structured to work out of the box
Once downloaded, all scripts should be able to be run as is, without any modification to the repository structure or scripts themselves, as long as the required packages and programs are installed.

#### Phyloseq objects used in the manuscript are included in this repository
As stated in the directory tree, phyloseq objects used in the manuscript for datasets 1 and 2 are located in the `PhyloseqObjects/` directory. Running scripts in `Dataset1_Scripts/`, `Dataset2_Scripts/`, and `Joint_Analyses_Scripts/` directories using these phyloseq objects should give same results as reported in the manuscript.

#### Results for analyses utilizing manuscript phyloseq objects are included in this repository
Results that are outputted by analyses scripts and that were reported in the manuscript are located in the `Script_Output/` directory. Any results that are not found in that directory should be in the supplementary material of the manuscript.

*WARNING: As the scripts are currently set up, the results that are included in the `Script_Output/` repository will be overwritten if re-running the analyses scripts without moving or renaming the result files first.*

#### There are two types of scripts stored in this repository
Scripts with the extension `.R` are R scripts written in R programming language used to perform differential abundance analyses and generate plots. Scripts that have a `.job` extension are shell scripts that were used to submit `.R` scripts to a SLURM scheduling system on a high performance computing cluster. Each `.R` script should have an accompanying `.job` script.

#### Each script should have an unfiltered and filtered data version
All differential abundance tests and concordance calculation scripts should have both an unfiltered data (no extra file suffix) and a filtered data (file suffix = "\_10%") version. This also goes for results found in the `Script_Output/` directory (except for hierarchical clustering and heatmap).

#### Perform your own comparison study
An R function `DAreport` has been added to the repository that allows comparisons detailed in the manuscript (Figures 1 and 3) to be performed on a user supplied feature abundance table. This was added in an attempt to address one of the limitations of the comparison study which is that DA methods will have varying performance depending on the underlying dataset being analyzed, therefore, `DAreport` was written in order to provide others with a user friendly way of running the comparison study on their own datasets to get comparison results specific to their dataset of interest. See the `README.md` file within the directory `DAreport/` for further details on how to run the function.
