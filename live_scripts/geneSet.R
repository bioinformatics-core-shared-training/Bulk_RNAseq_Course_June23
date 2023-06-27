library(tidyverse)
library(clusterProfiler)

search_kegg_organism('mouse', by = 'common_name')
# code is mmu

shrink.d11 <- readRDS("RObjects/Shrunk_Results.d11.rds")

sigGenes <- shrink.d11 %>%
  drop_na(Entrez, FDR) %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  pull(Entrez)

keggRes <- enrichKEGG(gene = sigGenes, organism = 'mmu')

as_tibble(keggRes)

browseKEGG(keggRes, 'mmu04612')

library(pathview)

logFC <- shrink.d11$logFC
names(logFC) <- shrink.d11$Entrez
pathview(gene.data = logFC,
         pathway.id = 'mmu04612',
         species = 'mmu',
         limit = list(gene=20, cpd=1))

# Exercise 1

logFC <- shrink.d11 %>%
  drop_na(FDR, Entrez) %>%
  filter(FDR < 0.01) %>%
  pull(logFC, Entrez)
  
pathview(gene.data = logFC,
         pathway.id = 'mmu04612',
         species = 'mmu',
         limit = list(gene=5, cpd=1))
  
# Go terms

library(org.Mm.eg.db)
sigGenes_GO <- shrink.d11 %>%
  drop_na(FDR) %>%
  filter(FDR < 0.01 & abs(logFC) > 2) %>%
  pull(GeneID)

universe <- shrink.d11$GeneID

ego <- enrichGO(gene = sigGenes_GO,
                universe = universe,
                OrgDb = org.Mm.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",
                pvalueCutoff = 0.01,
                readable = TRUE)
ego
as_tibble(ego)

barplot(ego, showCategory = 20)
dotplot(ego, font.size = 14)

library(enrichplot)

ego_pt <- pairwise_termsim(ego)
emapplot(ego_pt, cex_label_category = 0.25)

# GSEA
library(msigdbr)

rankedGenes <- shrink.d11 %>%
  filter(!is.na(GeneID)) %>%
  mutate(rank = logFC) %>%
  arrange(desc(rank)) %>%
  pull(rank, GeneID)

term2gene <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)
term2name <- msigdbr(species = "Mus musculus", category = "H") %>%
  dplyr::select(gs_name, gs_description) %>%
  distinct()

gseaRes <- GSEA(rankedGenes,
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pvalueCutoff = 1,
                minGSSize = 15,
                maxGSSize = 500)

as_tibble(gseaRes) %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, wt=-p.adjust) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digits=3)) %>%
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))

gseaplot(gseaRes,
         geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE",
         title = "HALLMARK_INFLAMMATORY_RESPONSE")

# Exercise 2

rankedGenes.e11 <- shrink.d11 %>%
  drop_na(GeneID, pvalue, logFC) %>%
  mutate(rank = -log10(pvalue) * sign(logFC)) %>%
  arrange(desc(rank)) %>%
  pull(rank, GeneID)

gseaRes.e11 <- GSEA(rankedGenes.e11,
                    TERM2GENE = term2gene,
                    TERM2NAME = term2name,
                    pvalueCutoff = 1,
                    minGSSize = 15,
                    maxGSSize = 500)

e11.tab <- as_tibble(gseaRes.e11) %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, wt=-p.adjust) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digits=3)) %>%
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))

gseaplot(gseaRes,
         geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
         title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION")

# day 33

shrink.d33 <- readRDS("RObjects/Shrunk_Results.d33.rds")

rankedGenes.e33 <- shrink.d33 %>%
  drop_na(GeneID, pvalue, logFC) %>%
  mutate(rank = -log10(pvalue) * sign(logFC)) %>%
  arrange(desc(rank)) %>%
  pull(rank, GeneID)

gseaRes.e33 <- GSEA(rankedGenes.e33,
                    TERM2GENE = term2gene,
                    TERM2NAME = term2name,
                    pvalueCutoff = 1,
                    minGSSize = 15,
                    maxGSSize = 500)

e33.tab <- as_tibble(gseaRes.e33) %>%
  arrange(desc(abs(NES))) %>%
  top_n(10, wt=-p.adjust) %>%
  mutate(across(c("enrichmentScore", "NES"), round, digits=3)) %>%
  mutate(across(c("pvalue", "p.adjust", "qvalue"), scales::scientific))

gseaplot(gseaRes,
         geneSetID = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
         title = "HALLMARK_INTERFERON_ALPHA_RESPONSE")

