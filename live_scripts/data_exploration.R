library(tximport)
library(DESeq2)
library(tidyverse)
# read in metadata
sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)
files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
files <- set_names(files, sampleinfo$SampleName)
tx2gene <- read_tsv("references/tx2gene.tsv")
# read in count data
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
saveRDS(txi, file="salmon_output/txi.rds")

tpm <- tximport(files, type="salmon",countsFromAbundance = "lengthScaledTPM",
                tx2gene=tx2gene)

rawCounts <- round(txi$counts, 0)
dim(rawCounts)

keep <- rowSums(rawCounts) >5
filtCounts <- rawCounts[keep,]
dim(filtCounts)

summary(filtCounts)

boxplot(filtCounts, main='Raw counts', las=2)
plot(rowMeans(filtCounts), rowSds(filtCounts),
     main='Raw coutns: sd vs mean', xlim=c(0,10000), ylim=c(0,5000))

# log2 transformation 
logcounts <- log2(filtCounts +1)

statusCols <- case_when(sampleinfo$Status=="Infected" ~ "red",
                        sampleinfo$Status=="Uninfected" ~ "orange")

boxplot(logcounts, xlab="", ylab="Log2(counts)", las=2,
        col=statusCols, main="Log2(Counts)")
abline(h=median(logcounts), col="blue")

plot(rowMeans(logcounts), rowSds(logcounts), main="Log2 Counts: sd vs mean")

# VST: variance stabilizing transformation
vst_counts <- vst(filtCounts)
boxplot(vst_counts, xlab="", ylab="VST counts", las=2, col=statusCols)
abline(h=median(vst_counts),col="blue")

plot(rowMeans(vst_counts), rowSds(vst_counts), main='VST counts: sd vs mean')

# PCA
library(ggfortify)

rlogcounts <- rlog(filtCounts)

pcDat <- prcomp(t(rlogcounts))
autoplot(pcDat)

autoplot(pcDat, data=sampleinfo, colour="Status", shape="TimePoint", size=5)
library(ggrepel)

autoplot(pcDat, data=sampleinfo, colour="Status", shape="TimePoint", size=5)+
  geom_text_repel(aes(x=PC1, y=PC2, label=SampleName),box.padding = 0.8)

sampleinfo <- mutate(sampleinfo, 
                     Status=case_when(SampleName=="SRR7657882" ~ "Uninfected",
                                      SampleName=="SRR7657873" ~ "Infected", 
                                      TRUE ~ Status))
write_tsv(sampleinfo, "results/SampleInfo_corrected.txt")

autoplot(pcDat, data=sampleinfo, colour="Status", shape="TimePoint", size=5)

# hierachical clustering
library(ggdendro)
hclDat <- t(rlogcounts) %>% 
  dist(method="euclidean") %>% 
  hclust()

ggdendrogram(hclDat, rotate=TRUE)

hclDat2 <- hclDat
hclDat2$labels <- str_c(sampleinfo$Status, ":", sampleinfo$TimePoint)
ggdendrogram(hclDat2, rotate=TRUE)


