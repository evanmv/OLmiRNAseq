library(readxl)
library(tidyverse)

miR129_miRDB <- read_xlsx("Data/miR-129-2-3p-miRDB.xlsx")
miR129_TS <- read.delim("Data/TargetScan8.0__miR-129-3p.Human.predicted_targets.txt")
merge <- filter(miR129_TS, miR129_TS$Target.gene %in% miR129_miRDB$`Gene Symbol`)
merge.f <- filter(miR129_TS, miR129_TS$Cumulative.weighted.context...score < -0.5 & miR129_TS$Aggregate.PCT > 0.5)
miR129_Targets <- merge.f$Target.gene

#for loop? or function?