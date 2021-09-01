###### Holly Chapapas - 1Sept2021 - TSCI5230 Lecture 2 Code #####

#Generate R report with Cntrl+Shift+K or click little notebook button

#+ init, echo=FALSE, message=FALSE, warning=FALSE,results='hide'
# debug <- 0;
# knitr::opts_chunk$set(echo=debug>0, warning=debug>0, message=debug>0);

# library(GGally)
# library(rio)
# library(dplyr)
# library(pander)
# library(synthpop)
# library(BiocManager)
# #BiocManager::install("recount")
# #BiocManager::install("GenomeInfoDbData")
# library(recount)
# BiocManager::install("microRNA")
# library(microRNA)
# BiocManager::install("miRNApath")
# library(miRNApath)

# #Libraries R currently sees
# search() %>% pander()

# Select dataset to work with
recount::download_study("SRP050223")
study <- load("SRP050223/rse_gene.Rdata")

# Examine the read count matrix
cts <- assay(rse_gene)

# Examine the sample metadata to see if it's ready for DESeq2 analysis
cd <- as.data.frame(colData(rse_gene))

# Use dat0 for purposes of class labeling this is the miRNA metadata
dat0 <- cd

#' Summarizes various facts
dat0cb <- codebook.syn(dat0)

# Create synthetic version of your data
dat0s <- syn(dat0)

#' Has all the same statistical properties as dat0 but nothing else
compare(dat0s,dat0)


########### Implement RNAseq workflow for next class until error is reached
