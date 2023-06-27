#!/usr/bin/Rscript
# Description: Day 3 of bulk RNA-seq Course: annotation and visualisation
# Load libraries
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)

# Annotation ----
## Load data ----
ddsObj.interaction     <- readRDS("RObjects/DESeqDataSet.interaction.rds")
results.interaction.11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
results.interaction.33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")

## Add the annotation to deseq2 results
### annotationhub
ah <- AnnotationHub()
ah
## Select the data that we need
MouseEnsDb <- query(ah, c('EnsDb', 'Mus musculus', '102'))[[1]]

## Extract genes
annotations <- genes(MouseEnsDb, return.type='data.frame')
annot <- annotations %>%
  dplyr::select(gene_id, gene_name, entrezid) %>%
  dplyr::filter(gene_id %in% rownames(results.interaction.11))

## query to an already made annotation file
ensemblAnnot <- readRDS("RObjects/Ensembl_annotations.rds")

## Annotate results ----
annot.interaction.11 <- as.data.frame(results.interaction.11)  %>%
  rownames_to_column('GeneID') %>%
  left_join(ensemblAnnot, 'GeneID') %>%
  rename(logFC=log2FoldChange, FDR=padj)

write_tsv(annot.interaction.11, "results/Interaction.11_Results_Annotated.txt")

# Visualisation ----
## sanity check
hist(annot.interaction.11$pvalue)

## Shrink values to avoid FP
ddsShrink.11 <- lfcShrink(ddsObj.interaction, 
                          res = results.interaction.11,
                          type = "ashr")

shrinkTab.11 <- as.data.frame(ddsShrink.11) %>%
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)

## MA plot ----
par(mfrow=c(1,2))
plotMA(results.interaction.11, alpha=0.05, main='raw')
plotMA(ddsShrink.11, alpha=0.05, main='shrinked')

## Volcano plot
volcano.Tab.11 <- shrinkTab.11 %>%
  mutate(`-log10(pvalue)`=-log10(pvalue))

ggplot(volcano.Tab.11, aes(x = logFC, y= `-log10(pvalue)`)) +
  geom_point(aes(col=FDR < 0.05), size=1)+ 
  geom_text_repel(data = ~top_n(.x, 4, wt=-FDR), aes(label=Symbol))


# Exercise 1 ----
## Shrink the results for the 33 days contrast
ddsShrink.33 <- lfcShrink(ddsObj.interaction, 
                          res = results.interaction.33,
                          type = "ashr")

## Create a new column of -log10(pvalue) values in your shrinkTab 
## for 33 days.
volcano.Tab.33 <- as.data.frame(ddsShrink.33) %>%
  mutate(`-log10(pvalue)`=-log10(pvalue))

## Create a plot with points coloured by FDR < 0.05 similar to 
## how we did in the first volcano plot
### Add annotation
shrinkTab.33 <- volcano.Tab.33 %>%
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, "GeneID") %>% 
  rename(logFC=log2FoldChange, FDR=padj)


ggplot(shrinkTab.33, aes(x = logFC, y= `-log10(pvalue)`)) +
  geom_point(aes(col=FDR < 0.05), size=1) + 
  geom_text_repel(data = ~top_n(.x, 4, wt=-FDR), aes(label=Symbol))

## Compare d11 vs d33
vol1 <- ggplot(volcano.Tab.11, aes(x = logFC, y= `-log10(pvalue)`)) +
  geom_point(aes(col=FDR < 0.05), size=1)+ 
  geom_text_repel(data = ~top_n(.x, 4, wt=-FDR), aes(label=Symbol)) + 
  ggtitle('Day 11')
vol2 <- ggplot(shrinkTab.33, aes(x = logFC, y= `-log10(pvalue)`)) +
  geom_point(aes(col=FDR < 0.05), size=1) + 
  geom_text_repel(data = ~top_n(.x, 4, wt=-FDR), aes(label=Symbol)) + 
  ggtitle('Day 33')
vol1 + vol2

# Exercise 2
# Recreate MAplot in ggplot
maTab.33 <- shrinkTab.33 %>%
  mutate(M = log2(baseMean))


ggplot(maTab.33, aes(x=M, y=logFC)) + 
  geom_point(aes(col=FDR<0.05), size =1) +
  scale_y_continuous(limits = c(-4,4), 
                     oob = scales::squish) + 
  theme_bw()


## Strip charts ----
# Extract ensembl id using the gene name
geneID <- dplyr::filter(shrinkTab.11, Symbol=='Il10ra') %>% 
  pull(GeneID)

plotCounts(dds = ddsObj.interaction, 
           gene = geneID, 
           intgroup = c('TimePoint', 'Status', 'Replicate'),
           returnData = T) %>%
  ggplot(aes(x=Status, y= log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) + 
  expand_limits(y=0) +
  labs(title='Normalised counts for Interleukin 10 receptor alpha (Il10ra)')


# Exercise 3 ----
## Stripchar with another gene, e.g. Jchain
geneID <- dplyr::filter(shrinkTab.11, Symbol=='Jchain') %>% 
  pull(GeneID)

plotCounts(dds = ddsObj.interaction, 
           gene = geneID, 
           intgroup = c('TimePoint', 'Status', 'Replicate'),
           returnData = T) %>%
  ggplot(aes(x=Status, y= log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) + 
  expand_limits(y=0) +
  labs(title='Normalised counts for Jchain')


# Venn Diagram ----

vennDat <- tibble(Geneid=rownames(results.interaction.11)) %>%
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 & 
           !is.na(results.interaction.11$padj) & 
           results.interaction.11$log2FoldChange > 0) %>% 
  
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 & 
           !is.na(results.interaction.11$padj) & 
           results.interaction.11$log2FoldChange < 0) %>%
  
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 & 
           !is.na(results.interaction.33$padj) & 
           results.interaction.33$log2FoldChange > 0) %>%
  
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 & 
           !is.na(results.interaction.33$padj) & 
           results.interaction.33$log2FoldChange < 0) 

## generate plot
ggvenn(vennDat, set_name_size = 4)


# Heatmap ----
sigGenes <- shrinkTab.11 %>%
  top_n(300, wt=-FDR) %>%
  pull('GeneID')

plotDat <- vst(ddsObj.interaction)[sigGenes,] %>%
  assay()
# get x scores
z.mat <- t(scale(t(plotDat), center = T, scale = T))
# change color
myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2,0,2), myPalette)
Heatmap(z.mat, name='z-score', 
        col = myRamp, 
        show_row_names = FALSE)

ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")], 
                        col = list(Status = c("Uninfected" = "darkgreen", 
                                              "Infected" = "palegreen"), 
                                   TimePoint = c("d11" = "lightblue", 
                                                 "d33" = "darkblue")))
Heatmap(z.mat, name='z-score', 
        col = myRamp, 
        show_row_names = FALSE, top_annotation = ha1)

saveRDS(annot.interaction.11, file="results/Annotated_Results.d11.rds")
saveRDS(shrinkTab.11, file="results/Shrunk_Results.d11.rds")
saveRDS(annot.interaction.33, file="results/Annotated_Results.d33.rds")
saveRDS(shrinkTab.33, file="results/Shrunk_Results.d33.rds")












