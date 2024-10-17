library(readxl)
library(tidyverse)
library(plyr)

#Function to read in data from miRDB and TargetScan, compare, and filter for score and conservation. Outputs tibble with gene names
#Usage mirNAME <- targetList('mirNAME in file')
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
miR1249 <- targetList('1249-3p')
miR1268 <- targetList('1268')
miR483 <- targetList('483-3p')
miR6511 <- targetList('6511-3p')
miR296 <- targetList('296-5p')
miR1306 <- targetList('1306-5p')
miR4835p <- targetList('483-5p')
miR486 <- targetList('486-5p')
miR3960 <- targetList('3960')
miR1247 <- targetList('1247-5p')

#Dataframe with all downregulated miRs and targets
df <- join_all(list(miR129, miR204, miR210, miR320, miR375, miR3065, miR1249, miR1268, miR6511, miR296, miR1306, miR483, miR4835p, miR486, miR3960, miR1247), by = 0, type = 'full')

write_csv(df, "Res/miR_targetList.csv")

df2 <- read_csv("Res/miR_targetList.csv")
