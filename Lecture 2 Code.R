###### Holly Chapapas - 1Sept2021 - TSCI5230 Lecture 2 Code #####

#Generate R report with Cntrl+Shift+K or click little notebook button

#+ init, echo=FALSE, message=FALSE, warning=FALSE,results='hide'
# debug <- 0;
# knitr::opts_chunk$set(echo=debug>0, warning=debug>0, message=debug>0);

library(GGally)
library(rio)
library(dplyr)
library(pander)
library(synthpop)
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


########### Implement DEseq workflow for next class until error is reached #############
# Add to the colData to give a column with experimental conditions
rse_gene$condition <- c(rep("CTRL", 109), rep("EXP", 291))

# Do DESeq2 analysis
# -- Make the dds
dds <- DESeqDataSet(rse_gene, design = ~condition)
# -- Get rlog --- use VST for large datasets
#rld <- rlog(dds)
vst <- vst(dds)
# -- plot PCA
plotPCA(vst)
# -- analyze
dds <- DESeq(dds)
# -- get results
res <- results(dds, contrast = c("condition", "CTRL", "EXP"))
# -- plotMA
plotMA(res)
# LFC shrink - only if necessary
resNorm <- lfcShrink(dds = dds, res = res, type = "normal", coef = 2)
# -- plotMA with resNorm
plotMA(resNorm)
# Make a DF
resdf <- as.data.frame(resNorm)
View(resdf)
# -- convert ENSG to gene symbol
ens2sym <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c("SYMBOL"))

# -- wrangle the resdf and join with ens2sym map
resdf <- resdf %>%
  rownames_to_column() %>%
  mutate(GENEID = gsub(rowname, pattern = "\\..+", replacement = "")) %>%
  dplyr::select(-rowname) %>%
  inner_join(y = ens2sym, by = "GENEID")
View(resdf)
# -- volcano plot
EnhancedVolcano(resdf, lab = resdf$SYMBOL, pCutoff = 1e-100,
                FCcutoff = 3,
                x = "log2FoldChange", y = "padj")
# -- Get the over-expressed genes
dat <- resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange > 2)
  
write_csv(dat, file = "over_expressed_genes.csv")

# Get the under-expressed genes
resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange < -2) %>%
  write_csv(file = "under_expressed_genes.csv")

#################### CREATE SYNTHETIC DATA ####################
#remove character columns for synthetic data creation
dat0 <- subset(dat, select = -c(GENEID,SYMBOL) )
#write new csv file for later use 
write_csv(dat0, file = "over_expressed_genes2.csv")
#create synthetic data set
sim0 <- syn(dat0)












######### more analyses below

## DIY over-representation analysis ##

# Get the over-expressed genes as a vector
over_expressed_genes <- resdf %>%
  dplyr::filter(padj < .01 & log2FoldChange > 2) %>%
  pull(SYMBOL)

# Get the gene sets and wrangle
gene_sets <- msigdbr(species = "Homo sapiens", category = "C5")
gene_sets <- gene_sets %>%
  dplyr::select(gs_name, gene_symbol)

# Run over-representation analysis
egmt <- enricher(gene = over_expressed_genes,
                 TERM2GENE = gene_sets)
edf <- as.data.frame(egmt)
View(edf)

# Plot results with clusterProfiler
dotplot(egmt)
barplot(egmt)

##########################
# using enricher() to find the "curated" gene sets ("C2")
## GSEA ##

# Discrete cutoff -- or continuous distribution?
EnhancedVolcano(resdf, lab = resdf$SYMBOL,
                pCutoff = 1e-2,
                FCcutoff = 2,
                x = "log2FoldChange", y = "padj")

# Adding a score for GSEA
resdf2 <- resdf %>%
  arrange(padj) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange))

# Deal with inf
resdf2 <- resdf2 %>%
  mutate(padj = case_when(padj == 0 ~ .Machine$double.xmin,
                          TRUE ~ padj)) %>%
  mutate(gsea_metric = -log10(padj) * sign(log2FoldChange))

# Remove NAs and order by GSEA
resdf2 <- resdf2 %>%
  filter(! is.na(gsea_metric)) %>%
  arrange(desc(gsea_metric))
View(resdf2)

# GSEA value histogram
hist(resdf2$gsea_metric, breaks = 100)

# Get the ranked GSEA vector
ranks <- resdf2 %>%
  select(SYMBOL, gsea_metric) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  deframe()

# Run GSEA
gseares <- GSEA(geneList = ranks, 
                TERM2GENE = gene_sets)
gsearesdf <- as.data.frame(gseares)
View(gsearesdf)

# Make summary plots
dotplot(gseares)

# Make GSEA plot for "GO_CELL_MATRIX_ADHESION"
gseaplot(gseares, geneSetID = "GO_CELL_MATRIX_ADHESION",
         title = "GO_CELL_MATRIX_ADHESION")

# Make GSEA plot for "GO_RIBOSOMAL_SUBUNIT"
gseaplot(gseares, geneSetID = "GO_RIBOSOMAL_SUBUNIT",
         title = "GO_RIBOSOMAL_SUBUNIT")

# Make GSEA plot for top and bottom results
# -- Get top 4 over-expressed pathways
top_pathways <- gsearesdf %>%
  top_n(n = 4, wt = NES) %>%
  pull(ID)
# -- Make gseaplot for each and return as list
top_pathway_plots <- lapply(top_pathways, function(pathway) {
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})
# -- Arrange with labels as a multi-panel plot
top_pathway_plot <- ggarrange(plotlist = top_pathway_plots,
                              ncol = 2, nrow = 2, labels = "AUTO")
# -- Save it
ggsave(top_pathway_plot, filename = "top_GSEA_up.png",
       height = 11, width = 18)
# -- Repeat steps with top 4 under-expressed pathways
bottom_pathways <- gsearesdf %>%
  top_n(n = 4, wt = -NES) %>%
  pull(ID)
bottom_pathway_plots <- lapply(bottom_pathways, function(pathway) {
  gseaplot(gseares, geneSetID = pathway, title = pathway)
})
bottom_pathway_plot <- ggarrange(plotlist = bottom_pathway_plots,
                                 ncol = 2, nrow = 2, labels = "AUTO")
ggsave(bottom_pathway_plot, filename = "top_GSEA_down.png",
       height = 11, width = 18)
