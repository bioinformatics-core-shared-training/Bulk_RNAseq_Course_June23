# Introduction to using DESeq2 Models 

# libraries ----
library(DESeq2)
library(tidyverse)


# Load the data and set up ----

# 1. the counts from salmon
txi <- readRDS('RObjects/txi.rds')
txi

# 2. sample meta-data

sampleinfo <- read_tsv('data/samplesheet_corrected.tsv', col_types = 'cccc')

all( sampleinfo$SampleName == colnames(txi$counts) )

# 3. set up the model 

sampleinfo


simple.model <- as.formula( ~ Status    )
simple.model

cbind(sampleinfo,   model.matrix(simple.model, data = sampleinfo) )


# because of alphabetical assignment our beta 0 / 1 are the wrong way around
sampleinfo <- mutate(sampleinfo, Status = fct_relevel(Status, "Uninfected"))

cbind(sampleinfo, model.matrix(simple.model, data = sampleinfo))



# Import the things into DESEq2 ----

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData  = sampleinfo,
                                       design = simple.model)
ddsObj.raw 


# subset the low count genes ----

keep <- rowSums(counts(ddsObj.raw)) > 5 
keep


ddsObj.filt <- ddsObj.raw[keep ,  ]
ddsObj.filt



# DESeq2 Manual ----

ddsObj <- ddsObj.filt


# step 1 - estimate size factors ----

ddsObj <- estimateSizeFactors(ddsObj)

# have a look at the original
normalizationFactors(ddsObj.filt)
# now for the normalised one
normalizationFactors(ddsObj)

colSums(normalizationFactors(ddsObj))

# have a look at the non-normalised counts
logCount <- log(counts(ddsObj, normalized = F   ) + 1 )
# plot an MA
limma::plotMA(logCount, array = 5, ylim = c(-5,5) )
abline(h = 0, col = 'red')


# now check for the normalised 
logNomalisedCounts <- log2(counts(ddsObj, normalized = T) +1)

limma::plotMA(logNomalisedCounts, array = 5, ylim = c(-5,5) )
abline(h = 0, col = 'red')


# Step 2 estimate dispertions  ----


ddsObj <- estimateDispersions(ddsObj)

plotDispEsts(ddsObj)



# Step 3 GLM and Wald Statistic ----

ddsObj <- nbinomWaldTest(ddsObj)



# DESEq 2 automatically ----


ddsObj.filt


# run all 3 steps 
ddsObj <- DESeq(ddsObj.filt)


# results function

results.simple <- results( ddsObj, alpha = 0.01   )
results.simple



# excercise 1 - Up and Down Regulated Genes ----

results.simple

sum(results.simple$padj < 0.05)



summary(results.simple)


# Try the additive model ----

# 1 . The gene counts from salmon 

txi <- readRDS('RObjects/txi.rds')

# 2. Meta - data for sample

sampleinfo <- read_tsv('data/samplesheet_corrected.tsv', col_types = 'cccc')
sampleinfo <- mutate(sampleinfo, Status = fct_relevel(Status, "Uninfected"))

# 3. The additive model 

additve.model <- as.formula(~ TimePoint + Status )

# ,make the deseq2 object

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additve.model)


keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep , ]

#Exercise 2 ----


ddsObj.filt <- DESeq(ddsObj.filt)


results.additive <- results(ddsObj.filt, alpha = 0.01)
results.additive

summary(results.additive)

resultsNames(ddsObj.filt)


# Default contrasts ----

cbind(sampleinfo, model.matrix(additve.model, data = sampleinfo))


results.InfectedVsUninfected <- results.additive
rm(results.additive)


# Exercise 3 ----


?results

results.d33Vsd11 <- results(ddsObj.filt,
                            name = 'TimePoint_d33_vs_d11',
                            alpha = 0.01)

summary(results.d33Vsd11)


# Quick look at the PCA
vstCount <- vst(ddsObj.filt, blind = T)
plotPCA(vstCount, intgroup = c('Status','TimePoint'))



# comparing models using the LRT ----

# comparing the simple vs additve ----

ddsO.LRT <- DESeq(ddsObj.filt, 
                  test = 'LRT',
                  reduced = simple.model)

results.Additve_v_simple <- results(ddsO.LRT)
results.Additve_v_simple 


sum(results.Additve_v_simple$padj < 0.01, na.rm = T)



# Exercise 4 - Comparing additive vs interaction ----

interaction.model <- as.formula( ~ TimePoint * Status   )

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi, 
                                       colData = sampleinfo,
                                       design = interaction.model)
keep <- rowSums(counts(ddsObj.raw)) >5
ddsObj.filt <- ddsObj.raw[keep,]


#run deseq
ddsObj.interaction <- DESeq(ddsObj.filt)



ddsO.LRT <- DESeq(ddsObj.interaction, 
                  test = 'LRT',
                  reduced = additve.model)



results.Interaction_vs_additive <- results(ddsO.LRT, alpha = 0.01)

table(results.Interaction_vs_additive$padj < 0.01)


resultsNames(ddsO.LRT)


results.interaction.11 <-results(ddsObj.interaction,
                                 name = "Status_Infected_vs_Uninfected",
                                 alpha = 0.05)

summary(results.interaction.11)


results.interaction.33 <- results(ddsObj.interaction,
                                  contrast = list(c(  "Status_Infected_vs_Uninfected",
                                                 "TimePointd33.StatusInfected"  )),
                                  alpha = 0.05)


summary(results.interaction.33)


# exercise 5 - Trying to get the DEGs for the timepoints ----



# 1 uninfected mice - d33 vs d11 

results.d33_v_d11_uninfected <- results(ddsObj.interaction,
                                        name = 'TimePoint_d33_vs_d11',
                                        alpha = 0.05)

summary(results.d33_v_d11_uninfected)


results.d33_v_d11_infected <- results(ddsObj.interaction,
                                      contrast = list(c('TimePoint_d33_vs_d11',
                                                        'TimePointd33.StatusInfected')))
summary(results.d33_v_d11_infected)
# 2 infected mice - d33 vs d11



#saving the tables ----

write_tsv(sampleinfo, 'results/SampleInfo_corrected.txt')
saveRDS(ddsObj.interaction, 'results/DESeqDataSet.interaction.rds')
saveRDS(results.interaction.11, 'results/DESeqResults.interaction_d11.rds')
saveRDS(results.interaction.33, 'results/DESeqResults.interaction_d33.rds')
















































































































































  















