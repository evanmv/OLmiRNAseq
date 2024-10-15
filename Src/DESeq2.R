##Packages -----
library(tidyverse)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(regionReport)
library(RColorBrewer)
library(gt)
library(gtExtras)

BiocManager::install("isomiRs")
library(isomiRs)
##Data -----
countMatrix <- read.delim("Data/GSE62809_miraligner_Count_matrix.txt")

#Set MicroRNA as row names and delete column
row.names(countMatrix) <- countMatrix$MicroRNA
countMatrix$MicroRNA <- NULL 

colData <- read.delim("Doc/studydesign.txt", row.names = 1)

#Rearrange columns to match study design (ascending order by sampleID) and remove MicroRNA
order <- sort(colnames(countMatrix))
countMatrix <- countMatrix[, order]

all(row.names(colData) == colnames(countMatrix))   

#Set NAs to 0
countMatrix[is.na(countMatrix)] <- 0
          
##DESeq2 -----
dds <- DESeqDataSetFromMatrix(countData = round(countMatrix),
                              colData = colData,
                              design = ~ group)

#Pre-filtering ??
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >=10) >= smallestGroupSize
dds <- dds[keep, ]

#DGE analysis
dds <- DESeq(dds)
res <- results(dds)
?results()

summary(res)
#Shrinkage of effect size with apeglm
resultsNames(dds)
resLFC <- lfcShrink(dds, 
                    coef = "group_progressive_vs_non.progressive",
                    type = "apeglm")
summary(resLFC)
head(resLFC)
plotMA(resLFC, ylim = c(-2,2))

#Order and filter for sig
resOrd <- res[order(res$pvalue),]
resSig <- subset(resOrd, pvalue < 0.05)
summary(resSig)
head(resSig)

sigMirIDs <- resSig@rownames
sigMirs <- as_tibble(resSig) %>%
  mutate(mirID = sigMirIDs, .before = 1)
sigMirsdown <-sigMirs %>%
  filter(sigMirs$log2FoldChange < -1)
sigMirsup <- sigMirs %>%
  filter(sigMirs$log2FoldChange > 1)
write_csv(sigMirs, "sigMirs.csv")  
write_csv(sigMirsup, "Res/sigMirsUp.csv")
write_csv(sigMirsdown, "Res/sigMirsDown.csv")

#Read in data tables of sig Mirs
sigMirs <- read.csv("Res/sigMirs.csv")
sigMirsDown <- read.csv("Res/sigMirsDown.csv")
sigMirsUp <- read.csv("Res/sigMirsUp.csv")

#GT table
#Heatmap vector of colors
myheatcolors <- brewer.pal(name="RdBu", n=11)
?brewer.pa1 ?filter
DEtable <- sigMirs %>% 
  mutate(Change = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>%
  filter(log2FoldChange >=1| log2FoldChange <= -1) %>%
  select(mirID, log2FoldChange, pvalue, Change) %>%
  gt(groupname_col = "Change") %>%
  tab_header(
    title = md("**DE miRNAs in progressive Leukoplakia v non-progressive**")
  ) %>%
  gt_color_rows(
    columns = log2FoldChange,
    domain = c(-3, 3),
    palette = rev(myheatcolors)
  ) %>%
  tab_row_group(
    label = md("**Up**"),
    rows = Change == "Up"
  ) %>%
  tab_row_group(
    label = md("**Down**"),
    rows = Change == "Down"
  ) %>%
  row_group_order(groups = c("**Up**", "**Down**"))

gtsave(DEtable, "Res/tableDEmirs2.png") #gtsave doesn't work on fox (no chrome)


##PCA -----
#Variance stabilizing transformation

vsd <- varianceStabilizingTransformation(dds)
head(assay(vsd), 3)

#PCA 
#plotPCA(vsd, intgroup="group") #use returnData=TRUE to save pca data to object
pcaData <- plotPCA(vsd, intgroup="group", returnData=TRUE)
percentVar <- formatC(100 * attr(pcaData, "percentVar"))

#plot w/ GGplot
PCA <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_point(size=5) +
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) +
  stat_ellipse(geom="polygon", level=0.95, alpha=0.2) + #Add elipses 
  coord_fixed() +
  theme_bw()

ggsave("Res/PCA2.png")
