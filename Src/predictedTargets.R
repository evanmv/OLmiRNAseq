library(readxl)
library(tidyverse)
install.packages("plyr")
library(plyr)

#Function to read in data from miRDB and TargetScan, compare, and filter for score and conservation. Outputs tibble with gene names
targetList <- function(Target) {
  miRDB <- read_xlsx(paste0("Data/miRDB/",Target,".xlsx"))
  
  miR_TS <- read.delim(paste0("Data/TargetScan/",Target,".txt"))
  
  merge <- miR_TS %>%
    filter(miR_TS$Target.gene %in% miRDB$'Gene Symbol')
  merge <- merge %>%
    filter(Cumulative.weighted.context...score < -0.5 & Aggregate.PCT > 0.5)
  
  Target <- merge$Target.gene %>%
    as_tibble() %>%
    `colnames<-`(Target)
    
}

miR129 <- targetList('129-2-3p')
miR204 <- targetList('204-5p')
miR210 <- targetList('210-5p')
miR320 <- targetList('320')
miR375 <- targetList('375')
miR3065 <- targetList('3065-5p')

x <- merge.data.frame(miR129, miR204, by = 0, all = T)[-1]

df <- join_all(list(miR129, miR204, miR210,miR320, miR375,miR3065), by = 0, type = 'full')

