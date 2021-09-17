## Overview

This directory contains an R function called `DAreport` to run various differential abundance (DA) methods compared in *Wallen 2021 BMC Bioinformatics* on a feature abundance table and generate a report showing the distribution of concordances between DA method results (Figure 1 of *Wallen 2021 BMC Bioinf*) and what differentially abundant features are being detected in relation to their mean relative abundance and fold change between two groups of interest (Figure 3 of *Wallen 2021 BMC Bioinf*). Results of DA methods and other data used in the reports will be returned in a list.

This function was created to address one of the main limitations of the comparison study, which is that DA methods will have varying performance depending on the underlying dataset being analyzed, therefore, `DAreport` was written in order to provide others with a user friendly way of running the comparison study on their own datasets to get comparison results specific to their dataset of interest.

An R program `run_DAreport.R` is also provided here to run `DAreport` on the command line, and will output the same reports as the function along with a multi-tab excel of contents that would normally be returned as a list if just using the `DAreport` function in R.

**Note:** certain modifications have been made since the original BMC Bioinformatics comparison study.
+ Robust CLR with matrix completion is now using default parameters for the `OptSpace` function.
+ ANCOM v2 R code has been updated to ANCOM v2.1 R code. See [FrederickHuangLin/ANCOM/scripts/ancom_v2.1.R](https://github.com/FrederickHuangLin/ANCOM/blob/master/scripts/ancom_v2.1.R). The code for ANCOM is contained within the function so it does not need to be sourced from the R script provided in the GitHub repository.
+ The strategy and functions used for GLM NBZI has changed. Now using `glm.nb` function from `stats` to perform negative binomial GLM and `zeroinfl` function from the `pscl` R package to perform zero-inflated negative binomial GLM. This helps with computation time for larger datasets. `NA` p-values that result from either function are given a 1.
+ Wilcoxon rank-sum test is now being used instead of Kruskal-Wallis as this would be the more appropriate test to use with only two groups of interest.

## Required packages

The following R packages are required to run `DAreport`. These packages do not need to be loaded prior to running `DAreport` as they will be loaded by the function internally.
```
phyloseq
dplyr
ROptSpace
ALDEx2
ANCOMBC
compositions
nlme
DESeq2
pscl
MASS
samr
edgeR
baySeq
metagenomeSeq
limma
ggplot2
reshape2
ggh4x
ggpubr
grid
scales
```
**Note:** make sure LEfSe scripts are accessible via your `PATH` variable before running `DAreport`. I found it easiest to download LEfSe by making its own conda environment and running `conda install -c biobakery lefse` to download it and its dependencies. Then other required packages can be installed via R using the R that was downloaded with LEfSe.

## Parameters

Both the R function and the stand-alone program have the same parameters that are required to run `DAreport`.

+ `ps` - Input phyloseq object that contains the observed, un-transformed or -normalized (to total sequence depth) abundances of features (taxa, gene families, pathways, etc.). Should have an `otu_table` component with the abundances, and a `sample_data` component with a variable defining the two groups of interest. This variable should be coded as `1` (for treatment/experimental group) and `0` (for control group). If running the `run_DAreport.R` program, this should be a phyloseq object saved as a `.rds` file.
+ `group.var` - The variable name in the `sample_data` component of the phyloseq object that contains data on which samples belong to which of the two groups of interest. Note only two group variables are accepted.
+ `skip` - Name(s) of DA method(s) that you would like to skip running. Must be one or more of `aldex2`, `ancom`, `ancombc`, `deseq2_wald`, `glm_nbzi`, `lefse`, `samseq`, `edger_exactTest_tmm`, `edger_exactTest_rle`, `bayseq`, `fitzig`, `fitfeat`, `limma_voom`, `t_tss`, `t_clr`, `t_rclr`, `w_tss`, `w_clr`, `w_rclr`, `glm_tss`, `glm_clr`, `glm_rclr`. If running the `DAreport` function in R, this should be a vector of method names. If running the `run_DAreport.R` program, method names should be supplied in a comma separated list with no spaces. Default is `NULL`, no skipping. Not required if running R function, but is when running the R program. Two methods you may consider skipping if analyzing large datasets are ANCOM and SAMseq.
+ `seed` - Numeric value to pass to the set.seed() function in order to to make reports consistent. Default is `1234`. Not required when running the R function, but is when running the R program.

## Returned data

+ `DAreport_concordance.pdf` - Visual report of the distributions of pairwise concordances for each DA method, the proportion of DA features detected by each method, and the relationship found between concordance and proportion of DA features. See Figure 1 of *Wallen 2021 BMC Bioinf* for an example of what this looks like.
+ `DAreport_MRA_vs_FC.pdf` - Visual report of what differentially abundant features are being detected in relation to their mean relative abundance and fold change between two groups of interest. See Figure 3 of *Wallen 2021 BMC Bioinf* for an example of what this looks like.
+ If using R function:
    + list of method results (under same names as those listed in the `skip` parameter)
    + `concordances` - `data.frame` of pairwise concordance between each DA method.
    + `prop_DA` - `data.frame` listing each method and it's corresponding proportion of DA features detected.
    + `mra_fc_enrichment_tests` - `data.frame` listing each method and it's Fisher test results that are showcased in the `DAreport_MRA_vs_FC.pdf`.
+ If using the R program:
    + `DAreport_results.xlsx` - Multi-tab excel that contains the same data returned as a list with the `DAreport` function.

## Implementation

Example usage of the R function `DAreport` using dataset 1 phyloseq object from *Wallen 2021 BMC Bioinf*.
```
> library(phyloseq)
> source('DAreport.R')

# read in phyloseq object
> ps <- readRDS('../PhyloseqObjects/Dataset1/phyloseq.rds')

# collapse amplicon sequence variants to Genus level
> ps <- tax_glom(ps, taxrank="Genus", NArm=FALSE)

# code case/control variable as 1 for case and 0 for control
# note: does not matter here whether its coded as numeric 1/0 or character 1/0, it will
# be converted to factor within DAreport regardless
> sample_data(ps)$case_control <- dplyr::recode(sample_data(ps)$case_control, Case=1L, Control=0L)

# run DAreport
> da.report <- DAreport(ps, "case_control")

```

Example usage of the R program `run_DAreport.R` using dataset 2 phyloseq object from *Wallen 2021 BMC Bioinf*.
+ First jump into R to format phyloseq object and write out to `.rds` file.
```
> library(phyloseq)
> ps <- readRDS('../PhyloseqObjects/Dataset2/phyloseq.rds')
> ps <- tax_glom(ps, taxrank="Genus", NArm=FALSE)
> sample_data(ps)$case_control <- dplyr::recode(sample_data(ps)$case_control, Case=1L, Control=0L)
> saveRDS(ps, 'ps.rds')
> quit()
```
+ Now run the R program on the outputted phyloseq object, skipping running of SAMseq.
```
$ ./run_DAreport.R ps.rds case_control samseq 1234

```
