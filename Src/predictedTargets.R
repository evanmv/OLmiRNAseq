library(readxl)
library(tidyverse)
install.packages("plyr")
library(plyr)

miR129_miRDB <- read_xlsx("Data/miR-129-2-3p-miRDB.xlsx")
miR129_TS <- read.delim("Data/TargetScan8.0__miR-129-3p.Human.predicted_targets.txt")
merge <- filter(miR129_TS, miR129_TS$Target.gene %in% miR129_miRDB$`Gene Symbol`)
merge.f <- filter(merge, merge$Cumulative.weighted.context...score < -0.5 & merge$Aggregate.PCT > 0.5)
miR129_Targets <- merge.f$Target.gene

#for loop? or function?
write_lines(miR129_Targets, "Res/miR129_Targets.txt") #May be better to combine all DEmiRs into table and write out
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

