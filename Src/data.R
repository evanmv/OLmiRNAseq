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
          
dds <- DESeqDataSetFromMatrix(countData = round(countMatrix),
                              colData = colData,
                              design = ~ group)
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

plotMA(resLFC, ylim = c(-2,2))


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
?brewer.pa1
DEtable <- sigMirs %>% 
  mutate(Change = ifelse(log2FoldChange > 0, 'Up', 'Down')) %>%
  select(sigMirIDs, log2FoldChange, pvalue, Change) %>%
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

gtsave(DEtable, "Res/tableDEmirs.png") #gtsave doesn't work on fox (no chrome)
