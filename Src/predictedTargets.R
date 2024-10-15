library(readxl)
library(tidyverse)

miR129_miRDB <- read_xlsx("Data/miR-129-2-3p-miRDB.xlsx")
miR129_TS <- read.delim("Data/TargetScan8.0__miR-129-3p.Human.predicted_targets.txt")
merge <- filter(miR129_TS, miR129_TS$Target.gene %in% miR129_miRDB$`Gene Symbol`)
merge.f <- filter(merge, merge$Cumulative.weighted.context...score < -0.5 & merge$Aggregate.PCT > 0.5)
miR129_Targets <- merge.f$Target.gene

#for loop? or function?
write_lines(miR129_Targets, "Res/miR129_Targets.txt") #May be better to combine all DEmiRs into table and write out
targetList <- function(miRDB, TargetScan) {
  miRDB <- read_xlsx(miRDB)
  
  miR_TS <- read.delim(TargetScan)
  
  merge <- miR_TS %>%
    filter(miR_TS$Target.gene %in% miR129_miRDB$'Gene Symbol')
  merge <- merge %>%
    filter(Cumulative.weighted.context...score < -0.5 & Aggregate.PCT > 0.5)
  
  lbl_Targets <- merge$Target.gene
    
}

miR1292 <- targetList("Data/miR-129-2-3p-miRDB.xlsx", "Data/TargetScan8.0__miR-129-3p.Human.predicted_targets.txt")

