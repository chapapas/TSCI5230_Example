#' ---
#' title: "TSCI5230 Class Examples on Proper R Reports"
#' author: 'Chapapas, Holly'
#' abstract:  "Proper documentation and coding to create clean and organized R Reports and inputting data"
#' description: 'Manuscript'
#' clean: false
#' self_contained: true
#' number_sections: false
#' keep_md: true
#' fig_caption: true
#' css: production.css
#' output:
#'   html_document:
#'    code_folding: hide
#'    toc: true
#'    toc_float: true
#' ---

# Generate R report with Cntrl+Shift+K or click little notebook button
#' This will be interpreted as text in the R Report
#' # This will be interpreted as a section name in R Report
#' ## This will be a section subheader in the R Report

#+ init, echo=FALSE, message=FALSE, warning=FALSE,results='hide'
inputdata <- c(dat0='sim_data.csv');
debug <- 0;
knitr::opts_chunk$set(echo=debug>0, warning=debug>0, message=debug>0);

#' Load libraries
library(GGally);
library(rio);
library(dplyr);
library(pander);
library(synthpop);
library(BiocManager)
# #BiocManager::install("recount")
# #BiocManager::install("GenomeInfoDbData")
library(recount)
library(DESeq2)
library(EnsDb.Hsapiens.v86)
library(tibble)
# #BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(readr)

#' Here are the libraries R currently sees:
search() %>% pander();

#' Select dataset to work with
recount::download_study("SRP050223")
study <- load("SRP050223/rse_gene.Rdata")

#' Examine the read count matrix
cts <- assay(rse_gene)

#' Examine the sample metadata to see if it's ready for DESeq2 analysis
cd <- as.data.frame(colData(rse_gene))

#' ## DESeq2 Workflow to Get Overexpressed Genes
#' Add to the colData to give a column with experimental conditions
rse_gene$condition <- c(rep("CTRL", 109), rep("EXP", 291))

#' Do DESeq2 analysis
#' Make the DESeq Data set (DDS)
dds <- DESeqDataSet(rse_gene, design = ~condition)
#' Get the vst
vst <- vst(dds)
#' Plot PCA
plotPCA(vst)
#' Analyze data
dds <- DESeq(dds)
#' Get the DESeq results
res <- results(dds, contrast = c("condition", "CTRL", "EXP"))
#' Plot MA
plotMA(res)
#' LFC shrink - only if necessary
resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)
#' Plot MA with resNorm
plotMA(resNorm)
#' Make a data frame
resdf <- as.data.frame(resNorm)
#' Convert ENSG to gene symbol
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c("SYMBOL"))

#' Wrangle the resdf and join with ens2sym map
resdf <- resdf %>%
  rownames_to_column() %>%
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) %>%
  dplyr::select(-rowname) %>%
  inner_join(y = ens2sym, by = "GENEID")

#' Create a volcano plot
EnhancedVolcano(resdf, lab = resdf$SYMBOL, pCutoff = 1e-100,
                FCcutoff = 3,
                x = "log2FoldChange", y = "padj")
#' Get the over-expressed genes
# Use dat0 for purposes of class labeling this is the miRNA metadata
dat0 <- resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange > 2)

write_csv(dat0, file = "over_expressed_genes.csv")

# Get the under-expressed genes
resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange < -2) %>%
  write_csv(file = "under_expressed_genes.csv")

#' ## Create Synthetic Data
#' Remove character columns for synthetic data creation
dat0s <- subset(dat0, select = -c(GENEID,SYMBOL) )
#' Write new csv file for later use 
write_csv(dat0s, file = "over_expressed_genes_syn.csv")
#' Create synthetic data set
sim0 <- syn(dat0s)
?subset
#' Summarizes various facts
dat0cb <- codebook.syn(dat0s)

#' Create synthetic version of your data
dat0s <- syn(dat0s)

#' Has all the same statistical properties as dat0 but nothing else
compare(dat0s,dat0)

#' Load data
dat0 <- import(inputdata['dat0'])

#' Make a scatterplot matrix
ggpairs(dat0);
#' Set all the two-value columns to be TRUE/FALSE
dat1 <- mutate(dat0, across(where( function(xx) length(unique(xx))==2), as.logical));
#' Now try the scatterplot matrix again
ggpairs(dat1);

#' ## Importing Synthetic Data
orig0 <- import('C:/Users/chapa/Documents/UTHSCSA/Course Work/Spring 2021/TSCI5230_Example/over_expressed_genes.csv')
orig0cb <- codebook.syn(orig0);
sim0 <- syn(orig0);
compare(sim0,orig0);
export(sim0$syn,"sim_data.csv");

inputdata <- c(dat0='over_expressed_genes_sim.csv');
if(file.exists('local.config.R')) source('local.config.R', local = TRUE, echo = FALSE);
dat0 <- import(inputdata['dat0']);

#' ## Filter
#' 
#' filter(dat0, row == "")
#' Using pander to format as a table
#filter(dat0, padj == "") %>% pander(split.tables=Inf,split.cells=Inf)

#' ## Slice
#' 
#' Slice takes a section of the data 
dplyr::slice(dat0, 1:5) %>% pander(split.tables=Inf,split.cells=Inf)
dplyr::slice_sample(dat0, n=5, replace = TRUE)

#' ## Arrange
#' 
#' Will sort your data in the order you specify
arrange(dat0, padj) %>% head()

#' ## Summarize
#' 
#' Aggregates data - can be mean, median, mode, ect.
#' Use Args to identify details of function for syntax
#' args(summarize)
#' 
#' This first line will do a mean of a specific column
summarize(dat0,mean=mean(stat))
#' This line will do the mean across all columns
dat0 %>% subset(dat0, select = -c(GENEID,SYMBOL) ) %>%
  summarize(dat0,across(.cols=everything(),mean))
#' Pipeline the means across column back into the data
dat0 %>% subset(dat0, select = -c(GENEID,SYMBOL) ) %>%
  summarize(across(.cols=everything(),mean))

#' ## Group By
#'
#' Undo the log fold change for grouping
2^dat0$log2FoldChange %>% round %>% hist(breaks=1000, xlim=c(0,20))
2^dat0$log2FoldChange %>% round %>% pmin(10) %>% table
#' Create a new column for the Fold Change
dat0$FoldChange <- 2^dat0$log2FoldChange %>% round %>% pmin(10)
#Group data by new column
group_by(dat0, FoldChange) %>% summarize(across(.cols=everything(),list(mean=mean,sd=sd)))


#' ## Plotting
#' 
#' Order the data by lowest (aka significant) padj value
dat0_arrange <- arrange(dat0, padj)
#' Create a data set with the top 20 significant genes for plotting
top20_dat0 <- dat0_arrange[-c(21:3468), ]

#' Plotting with ggplot2 package
ggplot(top20_dat0) +
  geom_point(aes(x = SYMBOL, y = log2FoldChange)) +
  xlab("Genes") +
  ylab("log2FoldChange") +
  ggtitle("Top 20 DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))

#+ facet_grid(condition~.)

# To rename a git file: git mv OLDNAME NEWNAME
# packageStatus()$inst %>% subset(Status=='upgrade') %>% rownames %>% gsub('\\.1','',.)


