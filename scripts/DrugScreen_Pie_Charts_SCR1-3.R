# ================================================================
# Title: DrugScreen_Pie_Charts_SCR1-3
# Description: Code to produce the final pie charts showing the fraction of compounds being hits for the primary screens
# (Fig. 2C)
# ================================================================

# ---- Load libraries ----
library(tidyverse)
library(plyr)
library(ggsci)

# ---- Prepare the color palettes ---- 
uchicago_colors <- pal_uchicago("default", alpha = 0.7)(5)
custom_order_colors <- uchicago_colors[c(3, 4, 2)] 
custom_pal <- pal_uchicago("default", alpha = 0.7)(9)
custom_pal_ordered <- custom_pal[c(1, 8, 4, 5, 6, 3, 9, 2)] 

# ---- Import StratoMiner output files ---- 
SCR1_AstroD_Pears <- read.delim("./data/SCR1PKEAstroDAggData_Pearson_WithControls.txt")
SCR1_AstroF_Pears <- read.delim("./data/SCR1PKEAstroFAggData_Pearson_WithControls.txt")
SCR2_AstroD_Pears <- read.delim("./data/SCR2TargetMAstroDAggData_Pearson_WithControls.txt")
SCR2_AstroF_Pears <- read.delim("./data/SCR2TargetMAstroFAggData_Pearson_WithControls.txt")
SCR3_AstroD_Pears <- read.delim("./data/Scr3FDAAstroDAggData_Pearson_WithControls.txt")
SCR3_AstroF_Pears <- read.delim("./data/Scr3FDAAstroFAggData_Pearson_WithControls.txt")

# ---- Import compounds lists ----
SCR1_comp_list <- read.csv("./data/SCR1_PKE_library/SCR1_picklist_PKE_150621.csv")
SCR2_comp_list <- read.csv("./data/SCR2_TargetMol_library/Analysis_pictures_from_Image_Xpress/Plate_annotation/SCR2_picklist_targetmol_allplates.csv", sep=";")
load("./data/SCR3_PrestWNumb_to_ChemName_only_used_compounds.rda")
SCR3_comp_list <- PrestWNumb_to_ChemName
rm(PrestWNumb_to_ChemName)

## SCR1 ##

# ---- Aggregate to have name lists ----
SCR1_comp_list <- SCR1_comp_list %>% ddply("Sample.Name", nrow)

# ---- Remove useless columns, add astrocyte batch and join tables ----
SCR1_AstroD_Pears <- SCR1_AstroD_Pears %>% filter(Category == "SAMPLE") %>%
  filter(!(ReagentID %in% c("AG957", "BMS-3", "BMS-5", "Imatinib", "JSH 23", "SP600125", "TPCA-1", "WP1066"))) %>%
  mutate(ID = str_c(Compound, Concentration, sep = "_")) %>%
  select(ID, DistanceScore, ReagentID, Concentration)
SCR1_AstroF_Pears <- SCR1_AstroF_Pears %>% filter(Category == "SAMPLE") %>%
  filter(!(ReagentID %in% c("AG957", "BMS-3", "BMS-5", "Imatinib", "JSH 23", "SP600125", "TPCA-1", "WP1066"))) %>%
  mutate(ID = str_c(Compound, Concentration, sep = "_")) %>%
  select(ID, DistanceScore)

# ---- Merge astroD and astroF tables SCR1 ----
SCR1_AllData <- merge(SCR1_AstroD_Pears, SCR1_AstroF_Pears, by.x = "ID", by.y = "ID", suffixes = c("_D","_F"))
rm(SCR1_AstroD_Pears, SCR1_AstroF_Pears)

# ---- Generate a data frame with the information if the compound is a hit for one batch or for both at one specific concentration ----
SCR1_AllData <- SCR1_AllData %>% filter(!(DistanceScore_D>0.7 & DistanceScore_F>0.7)) # Filter out the compounds that are not a hit for both astrocytes
SCR1_AllData <- SCR1_AllData %>% mutate(Hit_D = ifelse(DistanceScore_D <= 0.7, TRUE, FALSE), Hit_F = ifelse(DistanceScore_F <= 0.7, TRUE, FALSE), 
                                        Hit_both = (Hit_D & Hit_F), Hit_one = (Hit_D | Hit_F))
SCR1_AllData_Any_Cn <- aggregate(.~ReagentID, SCR1_AllData[,-c(1, 2, 4, 5)], any)

# ---- Add the information to SCR1_Comp_list to know for each compound if they are a hit both both batches, only one or none ----
Hit_batch <- rep("Not a hit", nrow(SCR1_comp_list))
for(i in 1:nrow(SCR1_comp_list)){
  for(j in 1:nrow(SCR1_AllData_Any_Cn)){
    if((SCR1_comp_list$Sample.Name[i]  == SCR1_AllData_Any_Cn$ReagentID[j]) & (SCR1_AllData_Any_Cn$Hit_both[j]  == TRUE)){
      Hit_batch[i] <- "Hit for 2 batches"}
    if((SCR1_comp_list$Sample.Name[i]  == SCR1_AllData_Any_Cn$ReagentID[j]) & (SCR1_AllData_Any_Cn$Hit_both[j]  == FALSE) & (SCR1_AllData_Any_Cn$Hit_one[j]  == TRUE)){
      Hit_batch[i] <- "Hit for 1 batch"}
  }}
SCR1_comp_list <- SCR1_comp_list %>% mutate(Hit_batch = as.factor(Hit_batch)) %>% select(Sample.Name, Hit_batch)
SCR1_summary <- SCR1_comp_list %>% ddply("Hit_batch", nrow)
SCR1_summary <- SCR1_summary %>% mutate(Fraction_hit = (V1/nrow(SCR1_comp_list))*100) %>%
  mutate(Fraction_hit = round(Fraction_hit, digits=1))

# ---- Make the pie chart - ggplot2 ----
SCR1_summary <- SCR1_summary %>%
  arrange(desc(Hit_batch)) %>%
  mutate(lab.ypos=cumsum(Fraction_hit)-0.5*Fraction_hit)
SCR1_PieChart <- SCR1_summary %>% 
  ggplot(aes(x="", y=Fraction_hit, fill=Hit_batch))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = V1), size=5, color = "black")+
  scale_fill_manual(values = custom_order_colors) +
  theme_void()+
  ggtitle("PKE library")
ggsave(file="Fig1C_SCR1_PieChart.pdf", SCR1_PieChart, path = "./results/", height = 3, width = 4)
rm(SCR1_AllData, SCR1_AllData_Any_Cn, SCR1_comp_list, SCR1_PieChart, SCR1_summary)

## SCR2 ##

# ---- Aggregate to have name lists ----
SCR2_comp_list <- SCR2_comp_list %>% mutate(Sample.Name = drugs.plate1) %>%
  ddply("Sample.Name", nrow)

# ---- Remove useless columns, add astrocyte batch and join tables ----
SCR2_AstroD_Pears <- SCR2_AstroD_Pears %>% filter(Category == "SAMPLE") %>%
  mutate(ID = str_c(Compound, Concentration, sep = "_")) %>%
  select(ID, DistanceScore, ReagentID, Concentration)
SCR2_AstroF_Pears <- SCR2_AstroF_Pears %>% filter(Category == "SAMPLE") %>%
  mutate(ID = str_c(Compound, Concentration, sep = "_")) %>%
  select(ID, DistanceScore)

# ---- Merge astroD and astroF tables ----
SCR2_AllData <- merge(SCR2_AstroD_Pears, SCR2_AstroF_Pears, by.x = "ID", by.y = "ID", suffixes = c("_D","_F"))
rm(SCR2_AstroD_Pears, SCR2_AstroF_Pears)

# ---- Generate a data frame with the information if the compound is a hit for one batch or for both at one specific concentration ----
SCR2_AllData <- SCR2_AllData %>% filter(!(DistanceScore_D>0.7 & DistanceScore_F>0.7)) # Filter out the compounds that are not a hit for both astrocytes
SCR2_AllData <- SCR2_AllData %>% mutate(Hit_D = ifelse(DistanceScore_D <= 0.7, TRUE, FALSE), Hit_F = ifelse(DistanceScore_F <= 0.7, TRUE, FALSE), 
                                        Hit_both = (Hit_D & Hit_F), Hit_one = (Hit_D | Hit_F))
SCR2_AllData_Any_Cn <- aggregate(.~ReagentID, SCR2_AllData[,-c(1, 2, 4, 5)], any)

# ---- Add the information to SCR2_Comp_list to know for each compound if they are a hit both both batches, only one or none ----
Hit_batch <- rep("Not a hit", nrow(SCR2_comp_list))
for(i in 1:nrow(SCR2_comp_list)){
  for(j in 1:nrow(SCR2_AllData_Any_Cn)){
    if((SCR2_comp_list$Sample.Name[i]  == SCR2_AllData_Any_Cn$ReagentID[j]) & (SCR2_AllData_Any_Cn$Hit_both[j]  == TRUE)){
      Hit_batch[i] <- "Hit for 2 batches"}
    if((SCR2_comp_list$Sample.Name[i]  == SCR2_AllData_Any_Cn$ReagentID[j]) & (SCR2_AllData_Any_Cn$Hit_both[j]  == FALSE) & (SCR2_AllData_Any_Cn$Hit_one[j]  == TRUE)){
      Hit_batch[i] <- "Hit for 1 batch"}
  }}
SCR2_comp_list <- SCR2_comp_list %>% mutate(Hit_batch = as.factor(Hit_batch)) %>% select(Sample.Name, Hit_batch)
SCR2_summary <- SCR2_comp_list %>% ddply("Hit_batch", nrow)
SCR2_summary <- SCR2_summary %>% mutate(Fraction_hit = (V1/nrow(SCR2_comp_list))*100) %>%
  mutate(Fraction_hit = round(Fraction_hit, digits=1))

# ---- Make the pie chart - ggplot2 ----
SCR2_summary <- SCR2_summary %>%
  arrange(desc(Hit_batch)) %>%
  mutate(lab.ypos=cumsum(Fraction_hit)-0.5*Fraction_hit)

SCR2_PieChart <- SCR2_summary %>% 
  ggplot(aes(x="", y=Fraction_hit, fill=Hit_batch))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = V1), size=5, color = "black")+
  scale_fill_manual(values = custom_order_colors) +
  theme_void()+
  ggtitle("Anti-cancer drugs library")
ggsave(file="Fig1C_SCR2_PieChart.pdf", SCR2_PieChart, path = "./results/", height = 3, width = 4)
rm(SCR2_AllData, SCR2_AllData_Any_Cn, SCR2_comp_list, SCR2_PieChart, SCR2_summary)

## SCR3 ##

# ---- Aggregate to have name lists ----
SCR3_comp_list <- SCR3_comp_list %>% mutate(Sample.Name = chemical_name) %>% select(Sample.Name, therapeutic_class, therapeutic_effect, target_type, target_name, target_mechanism)

# ---- Remove useless columns, add astrocyte batch and join tables ----
SCR3_AstroD_Pears <- SCR3_AstroD_Pears %>% separate_wider_delim(XPlate, "_", names = c("Astrocyte", NA, "Concentration", NA)) %>% 
  filter(Category == "SAMPLE") %>%
  filter(!(ReagentID %in% c("Tyrphostin-9", "PKC-412", "Erbstatin analog", "Rottlerin", "Ro 31-8220 (mesylate)"))) %>%
  mutate(ID = str_c(Compound, Concentration, sep = "_")) %>%
  select(ID, DistanceScore, ReagentID, Concentration)
SCR3_AstroF_Pears <- SCR3_AstroF_Pears %>% separate_wider_delim(XPlate, "_", names = c("Astrocyte", NA, "Concentration", NA)) %>% 
  filter(Category == "SAMPLE") %>%
  filter(!(ReagentID %in% c("Tyrphostin-9", "PKC-412", "Erbstatin analog", "Rottlerin", "Ro 31-8220 (mesylate)"))) %>%
  mutate(ID = str_c(Compound, Concentration, sep = "_")) %>%
  select(ID, DistanceScore)

# ---- Merge astroD and astroF tables SCR3 ----
SCR3_AllData <- merge(SCR3_AstroD_Pears, SCR3_AstroF_Pears, by.x = "ID", by.y = "ID", suffixes = c("_D","_F"))
rm(SCR3_AstroD_Pears, SCR3_AstroF_Pears)

# ---- Generate a data frame with the information if the compound is a hit for one batch or for both at one specific concentration ----
SCR3_AllData <- SCR3_AllData %>% filter(!(DistanceScore_D>0.7 & DistanceScore_F>0.7)) # Filter out the compounds that are not a hit for both astrocytes
SCR3_AllData <- SCR3_AllData %>% mutate(Hit_D = ifelse(DistanceScore_D <= 0.7, TRUE, FALSE), Hit_F = ifelse(DistanceScore_F <= 0.7, TRUE, FALSE), 
                                        Hit_both = (Hit_D & Hit_F), Hit_one = (Hit_D | Hit_F))
SCR3_AllData_Any_Cn <- aggregate(.~ReagentID, SCR3_AllData[,-c(1, 2, 4, 5)], any)

# ---- Add the information to SCR3_Comp_list to know for each compound if they are a hit both both batches, only one or none ----
Hit_batch <- rep("Not a hit", nrow(SCR3_comp_list))
for(i in 1:nrow(SCR3_comp_list)){
  for(j in 1:nrow(SCR3_AllData_Any_Cn)){
    if((SCR3_comp_list$Sample.Name[i]  == SCR3_AllData_Any_Cn$ReagentID[j]) & (SCR3_AllData_Any_Cn$Hit_both[j]  == TRUE)){
      Hit_batch[i] <- "Hit for 2 batches"}
    if((SCR3_comp_list$Sample.Name[i]  == SCR3_AllData_Any_Cn$ReagentID[j]) & (SCR3_AllData_Any_Cn$Hit_both[j]  == FALSE) & (SCR3_AllData_Any_Cn$Hit_one[j]  == TRUE)){
      Hit_batch[i] <- "Hit for 1 batch"}
  }}
SCR3_comp_list <- SCR3_comp_list %>% mutate(Hit_batch = as.factor(Hit_batch)) 
SCR3_summary <- SCR3_comp_list %>% select(Sample.Name, Hit_batch) %>% ddply("Hit_batch", nrow)
SCR3_summary <- SCR3_summary %>% mutate(Fraction_hit = (V1/nrow(SCR3_comp_list))*100) %>%
  mutate(Fraction_hit = round(Fraction_hit, digits=1))

# ---- Make the pie chart - ggplot2 ----
SCR3_summary <- SCR3_summary %>%
  arrange(desc(Hit_batch)) %>%
  mutate(lab.ypos=cumsum(Fraction_hit)-0.5*Fraction_hit)

SCR3_PieChart <- SCR3_summary %>% 
  ggplot(aes(x="", y=Fraction_hit, fill=Hit_batch))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = V1), size=5, color = "black")+
  scale_fill_manual(values = custom_order_colors) +
  theme_void()+
  ggtitle("FDA-approved drugs library")
ggsave(file="Fig1C_SCR3_PieChart.pdf", SCR3_PieChart, path = "./results/", height = 3, width = 4)

# ---- Make pie charts for SCR3 hits - therapeutic classes ----
SCR3_hits_only <- SCR3_comp_list %>% filter(Hit_batch == "Hit for 2 batches")

# ---- Add therapeutic categories with those <4 classed as "others" ----
SCR3_hits_only <- SCR3_hits_only %>% mutate(Therapeutic_group = NA)
for(i in 1:nrow(SCR3_hits_only)){
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Dermatology", "Diagnostic", "Gastroenterology", "Hematology", "Neuromuscular", "Ophthalmology", "Respiratory")){
    SCR3_hits_only$Therapeutic_group[i] <- "Others"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Allergology")){
    SCR3_hits_only$Therapeutic_group[i] <- "Allergology"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Cardiovascular")){
    SCR3_hits_only$Therapeutic_group[i] <- "Cardiovascular"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Central Nervous System")){
    SCR3_hits_only$Therapeutic_group[i] <- "Central Nervous System"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Endocrinology")){
    SCR3_hits_only$Therapeutic_group[i] <- "Endocrinology"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Infectiology")){
    SCR3_hits_only$Therapeutic_group[i] <- "Infectiology"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Metabolism")){
    SCR3_hits_only$Therapeutic_group[i] <- "Metabolism"}
  if(SCR3_hits_only$therapeutic_class[i] %in% c("Oncology")){
    SCR3_hits_only$Therapeutic_group[i] <- "Oncology"}
}

SCR3_summary_Therapeutic_group <- SCR3_hits_only %>% ddply("Therapeutic_group", nrow)
SCR3_summary_Therapeutic_group <- SCR3_summary_Therapeutic_group %>% mutate(Fraction_class = (V1/nrow(SCR3_hits_only))*100) %>%
  mutate(Fraction_class = round(Fraction_class, digits=1))
SCR3_summary_Therapeutic_group <- SCR3_summary_Therapeutic_group %>%
  arrange(desc(Therapeutic_group)) %>%
  mutate(lab.ypos=cumsum(Fraction_class)-0.5*Fraction_class)

SCR3_PieChart_ther_group <- SCR3_summary_Therapeutic_group %>% 
  ggplot(aes(x="", y=Fraction_class, fill=Therapeutic_group))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal_ordered) +
  theme_void()+
  ggtitle("Therapeutic classes")
ggsave(file="Fig1C_SCR3_PieChart_ther_class.pdf", SCR3_PieChart_ther_group, path = "./results/", height = 3, width = 4)



