
library(dplyr)

setwd("~/Dropbox/My documents/Projects/Future Fellowship/RESEARCH/NESP-coral-analysis-revised/SuppTable_CoralMapping")

coral_cat <- read.csv("coral_species_grouping ET_Genus GF.csv")
coral_desc <- read.csv("coral_data_mrt_GenusGF ET 9Nov22.csv")

SuppTab <- coral_cat %>% right_join(coral_desc) %>% 
  dplyr::select(Genus_GF, original_label = label, Species) %>%
  group_by(Genus_GF, original_label, Species) %>%
  summarize()

write.csv(SuppTab, "SuppTable.csv", row.names = F)

# Add Species_Comment to Table 2

tab2 <- read.csv("Table2_Indicator taxa_v2.csv")

tab2 <- tab2 %>% left_join(coral_desc, by = c("Indicator.taxa" = "Genus_GF"))

write.csv(tab2, "Tab2_withSpecies.csv", row.names = F)
