# ================================================================
# Title: DrugScreen_Hit_List
# Description: Code to produce the final hit lists and the pie chart to visualize their distribution
# (Fig 2F, Supp. 5B, Supp. table 6)
# ================================================================

# ---- Load libraries ----
library(tidyverse)
library(plyr)
library(ggsci)

# ---- Import StratoMiner output files ---- 
AstroG_Pears <- read.delim("./data/SCR4ConfScreenAstroG_allHits.txt")
AstroH_Pears <- read.delim("./data/SCR4ConfScreenAstroH_allHits.txt")
AstroI_Pears <- read.delim("./data/SCR4ConfScreenAstroI_allHits.txt")
AstroJ_Pears <- read.delim("./data/SCR4ConfScreenAstroJ_allHits.txt")

# ---- Prepare the color palettes ---- 
custom_pal <- pal_uchicago("default", alpha = 0.7)(9)
custom_pal_ordered <- custom_pal[c(1, 4, 2, 3)] 
custom_pal_ordered_2 <- custom_pal[c(1, 4)] 
custom_pal_ordered_3 <- custom_pal[c(4, 5, 2)] 
custom_pal_ordered_4 <- custom_pal[c(4, 6, 3, 9, 2)]
custom_pal_2 <- pal_uchicago("default", alpha = 0.9)(9)

# ---- Clean the data ---- 
wells <- c("n14", "n15", "o14", "o15", "b16", "b17", "c16", "c17")
AstroG_Pears <- AstroG_Pears %>% filter(!(Compound == "Procainamide (hydrochloride)" & WellLocation %in% wells))
AstroH_Pears <- AstroH_Pears %>% filter(!(Compound == "Procainamide (hydrochloride)" & WellLocation %in% wells))
AstroI_Pears <- AstroI_Pears %>% filter(!(Compound == "Procainamide (hydrochloride)" & WellLocation %in% wells))
AstroJ_Pears <- AstroJ_Pears %>% filter(!(Compound == "Procainamide (hydrochloride)" & WellLocation %in% wells))
AstroG_Pears <- AstroG_Pears %>% select(Compound, Concentration, Solvent, Category, IRstatut, XPlate, Well, DistanceScore,  pvalue) %>% 
  mutate(Astrocyte = "AstroG", Compound_Cn = str_c(Compound, Concentration, sep = "_") )
AstroH_Pears <- AstroH_Pears %>% select(Compound, Concentration, Solvent, Category, IRstatut, XPlate, Well, DistanceScore,  pvalue) %>% 
  mutate(Astrocyte = "AstroH", Compound_Cn = str_c(Compound, Concentration, sep = "_") )
AstroI_Pears <- AstroI_Pears %>% select(Compound, Concentration, Solvent, Category, IRstatut, XPlate, Well, DistanceScore,  pvalue) %>% 
  mutate(Astrocyte = "AstroI", Compound_Cn = str_c(Compound, Concentration, sep = "_") )
AstroJ_Pears <- AstroJ_Pears %>% select(Compound, Concentration, Solvent, Category, IRstatut, XPlate, Well, DistanceScore,  pvalue) %>% 
  mutate(Astrocyte = "AstroJ", Compound_Cn = str_c(Compound, Concentration, sep = "_") )

# ---- Separate in 2 df for each frames for the duplicate ---- 
AstroG_A <- AstroG_Pears %>% filter(str_detect(AstroG_Pears$Well,"[A-Z]2$|[A-Z]4$|[A-Z]6$|[A-Z]8$|[A-Z]10$|[A-Z]12$|[A-Z]14$|[A-Z]16$|[A-Z]18$|[A-Z]20$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))
AstroG_B <- AstroG_Pears %>% filter(str_detect(AstroG_Pears$Well,"[A-Z]3$|[A-Z]5$|[A-Z]7$|[A-Z]9$|[A-Z]11$|[A-Z]13$|[A-Z]15$|[A-Z]17$|[A-Z]19$|[A-Z]21$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))

AstroH_A <- AstroH_Pears %>% filter(str_detect(AstroH_Pears$Well,"[A-Z]2$|[A-Z]4$|[A-Z]6$|[A-Z]8$|[A-Z]10$|[A-Z]12$|[A-Z]14$|[A-Z]16$|[A-Z]18$|[A-Z]20$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))
AstroH_B <- AstroH_Pears %>% filter(str_detect(AstroH_Pears$Well,"[A-Z]3$|[A-Z]5$|[A-Z]7$|[A-Z]9$|[A-Z]11$|[A-Z]13$|[A-Z]15$|[A-Z]17$|[A-Z]19$|[A-Z]21$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))

AstroI_A <- AstroI_Pears %>% filter(str_detect(AstroI_Pears$Well,"[A-Z]2$|[A-Z]4$|[A-Z]6$|[A-Z]8$|[A-Z]10$|[A-Z]12$|[A-Z]14$|[A-Z]16$|[A-Z]18$|[A-Z]20$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))
AstroI_B <- AstroI_Pears %>% filter(str_detect(AstroI_Pears$Well,"[A-Z]3$|[A-Z]5$|[A-Z]7$|[A-Z]9$|[A-Z]11$|[A-Z]13$|[A-Z]15$|[A-Z]17$|[A-Z]19$|[A-Z]21$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))

AstroJ_A <- AstroJ_Pears %>% filter(str_detect(AstroJ_Pears$Well,"[A-Z]2$|[A-Z]4$|[A-Z]6$|[A-Z]8$|[A-Z]10$|[A-Z]12$|[A-Z]14$|[A-Z]16$|[A-Z]18$|[A-Z]20$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))
AstroJ_B <- AstroJ_Pears %>% filter(str_detect(AstroJ_Pears$Well,"[A-Z]3$|[A-Z]5$|[A-Z]7$|[A-Z]9$|[A-Z]11$|[A-Z]13$|[A-Z]15$|[A-Z]17$|[A-Z]19$|[A-Z]21$")) %>%
  filter(!(Compound %in% c("Water", "Medium DMSO", "High DMSO")))

# ---- Create a new df by merging the rep ---- 
AstroG_B <- AstroG_B %>% select(Compound_Cn, Compound, Concentration, IRstatut, Well, DistanceScore, pvalue, Astrocyte)
AstroG_rep <- merge(AstroG_A, AstroG_B, by = "Compound_Cn", suffixes = c("_rep1","_rep2"), all = TRUE)
 for (i in 1:nrow(AstroG_rep)){
   if(is.na(AstroG_rep$Compound_rep1[i])){
     AstroG_rep$Compound_rep1[i] <- AstroG_rep$Compound_rep2[i] 
     AstroG_rep$Concentration_rep1[i] <- AstroG_rep$Concentration_rep2[i]
     AstroG_rep$IRstatut_rep1[i] <- AstroG_rep$IRstatut_rep2[i]
     AstroG_rep$Astrocyte_rep1[i] <- AstroG_rep$Astrocyte_rep2[i]}}
AstroG_rep <- AstroG_rep %>% mutate(Compound= Compound_rep1, Concentration=Concentration_rep1, IR_statut=IRstatut_rep1, Astrocyte=Astrocyte_rep1) %>%
  select(Compound_Cn, Compound, Concentration, IR_statut, Astrocyte, DistanceScore_rep1, DistanceScore_rep2)

AstroH_B <- AstroH_B %>% select(Compound_Cn, Compound, Concentration, IRstatut, Well, DistanceScore, pvalue, Astrocyte)
AstroH_rep <- merge(AstroH_A, AstroH_B, by = "Compound_Cn", suffixes = c("_rep1","_rep2"), all = TRUE)
for (i in 1:nrow(AstroH_rep)){
  if(is.na(AstroH_rep$Compound_rep1[i])){
    AstroH_rep$Compound_rep1[i] <- AstroH_rep$Compound_rep2[i] 
    AstroH_rep$Concentration_rep1[i] <- AstroH_rep$Concentration_rep2[i]
    AstroH_rep$IRstatut_rep1[i] <- AstroH_rep$IRstatut_rep2[i]
    AstroH_rep$Astrocyte_rep1[i] <- AstroH_rep$Astrocyte_rep2[i]}}
AstroH_rep <- AstroH_rep %>% mutate(Compound= Compound_rep1, Concentration=Concentration_rep1, IR_statut=IRstatut_rep1, Astrocyte=Astrocyte_rep1) %>%
  select(Compound_Cn, Compound, Concentration, IR_statut, Astrocyte, DistanceScore_rep1, DistanceScore_rep2)

AstroI_B <- AstroI_B %>% select(Compound_Cn, Compound, Concentration, IRstatut, Well, DistanceScore, pvalue, Astrocyte)
AstroI_rep <- merge(AstroI_A, AstroI_B, by = "Compound_Cn", suffixes = c("_rep1","_rep2"), all = TRUE)
for (i in 1:nrow(AstroI_rep)){
  if(is.na(AstroI_rep$Compound_rep1[i])){
    AstroI_rep$Compound_rep1[i] <- AstroI_rep$Compound_rep2[i] 
    AstroI_rep$Concentration_rep1[i] <- AstroI_rep$Concentration_rep2[i]
    AstroI_rep$IRstatut_rep1[i] <- AstroI_rep$IRstatut_rep2[i]
    AstroI_rep$Astrocyte_rep1[i] <- AstroI_rep$Astrocyte_rep2[i]}}
AstroI_rep <- AstroI_rep %>% mutate(Compound= Compound_rep1, Concentration=Concentration_rep1, IR_statut=IRstatut_rep1, Astrocyte=Astrocyte_rep1) %>%
  select(Compound_Cn, Compound, Concentration, IR_statut, Astrocyte, DistanceScore_rep1, DistanceScore_rep2)

AstroJ_B <- AstroJ_B %>% select(Compound_Cn, Compound, Concentration, IRstatut, Well, DistanceScore, pvalue, Astrocyte)
AstroJ_rep <- merge(AstroJ_A, AstroJ_B, by = "Compound_Cn", suffixes = c("_rep1","_rep2"), all = TRUE)
for (i in 1:nrow(AstroJ_rep)){
  if(is.na(AstroJ_rep$Compound_rep1[i])){
    AstroJ_rep$Compound_rep1[i] <- AstroJ_rep$Compound_rep2[i] 
    AstroJ_rep$Concentration_rep1[i] <- AstroJ_rep$Concentration_rep2[i]
    AstroJ_rep$IRstatut_rep1[i] <- AstroJ_rep$IRstatut_rep2[i]
    AstroJ_rep$Astrocyte_rep1[i] <- AstroJ_rep$Astrocyte_rep2[i]}}
AstroJ_rep <- AstroJ_rep %>% mutate(Compound= Compound_rep1, Concentration=Concentration_rep1, IR_statut=IRstatut_rep1, Astrocyte=Astrocyte_rep1) %>%
  select(Compound_Cn, Compound, Concentration, IR_statut, Astrocyte, DistanceScore_rep1, DistanceScore_rep2)
rm(AstroG_A, AstroG_B, AstroH_A, AstroH_B, AstroI_A, AstroI_B, AstroJ_A, AstroJ_B, AstroG_Pears, AstroH_Pears, AstroI_Pears, AstroJ_Pears)

# ---- Calculate the mean Distance score ---- 
AstroG_rep <- AstroG_rep %>% rowwise() %>% mutate(Mean_DistanceScore = mean(c(DistanceScore_rep1, DistanceScore_rep2), na.rm = TRUE))
AstroH_rep <- AstroH_rep %>% rowwise() %>% mutate(Mean_DistanceScore = mean(c(DistanceScore_rep1, DistanceScore_rep2), na.rm = TRUE))
AstroI_rep <- AstroI_rep %>% rowwise() %>% mutate(Mean_DistanceScore = mean(c(DistanceScore_rep1, DistanceScore_rep2), na.rm = TRUE))
AstroJ_rep <- AstroJ_rep %>% rowwise() %>% mutate(Mean_DistanceScore = mean(c(DistanceScore_rep1, DistanceScore_rep2), na.rm = TRUE))

# ---- Merge the df for all batches and calculate the mean distance score ---- 
preAstroG <- AstroG_rep %>% mutate(Mean_DistanceScore_G = Mean_DistanceScore) %>% 
  select(Compound_Cn, Compound, Concentration, Astrocyte, IR_statut, Mean_DistanceScore_G) 
preAstroH <- AstroH_rep %>% mutate(Mean_DistanceScore_H = Mean_DistanceScore) %>% 
  select(Compound_Cn, Mean_DistanceScore_H)
preAstroI <- AstroI_rep %>% mutate(Mean_DistanceScore_I = Mean_DistanceScore) %>% 
  select(Compound_Cn, Mean_DistanceScore_I)
preAstroJ <- AstroJ_rep %>% mutate(Mean_DistanceScore_J = Mean_DistanceScore) %>% 
  select(Compound_Cn, Mean_DistanceScore_J)
pre_All_batches_1 <- merge(preAstroG, preAstroH, by = "Compound_Cn", all = TRUE)
pre_All_batches_2 <- merge(pre_All_batches_1, preAstroI, by = "Compound_Cn", all = TRUE)
All_batches <- merge(pre_All_batches_2, preAstroJ, by = "Compound_Cn", all = TRUE)
rm(preAstroG, preAstroH, preAstroI, preAstroJ, pre_All_batches_1, pre_All_batches_2)

# ---- Fix the compound names and concentration with no values ----
All_batches <- All_batches %>% select(-Compound, -Concentration) 
All_batches[c("Compound", "Concentration")] <- str_split_fixed(All_batches$Compound_Cn, "_", 2)

# ---- Calculate Mean Distance Score across batches ---- 
All_batches <- All_batches %>% rowwise() %>% mutate(Mean_DistanceScore_AllB = mean(c(Mean_DistanceScore_G, Mean_DistanceScore_H, Mean_DistanceScore_I, Mean_DistanceScore_J), na.rm = TRUE)) 
All_batches <- All_batches %>% select(Compound_Cn, Compound, Concentration, IR_statut, Mean_DistanceScore_AllB, Mean_DistanceScore_G, Mean_DistanceScore_H, Mean_DistanceScore_I, Mean_DistanceScore_J)

# ---- Calculate Mean Distance Score without AstroJ ---- 
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
} # Detach plyr if loaded as it inteferes with following script line
All_batches <- All_batches %>% rowwise() %>% mutate(Mean_DistanceScore_AllB_no_J = mean(c(Mean_DistanceScore_G, Mean_DistanceScore_H, Mean_DistanceScore_I), na.rm = TRUE))

# ---- Add a new column counting the number of NA values for each row ---- 
All_batches <- All_batches %>%
  rowwise() %>%
  mutate(NA_Count = sum(is.na(c(Mean_DistanceScore_G, 
                                Mean_DistanceScore_H, 
                                Mean_DistanceScore_I, 
                                Mean_DistanceScore_J)))) %>%
  ungroup()
# Reattach plyr
if (!"plyr" %in% loadedNamespaces()) {
  library(plyr)
}

# ---- Filter out the compounds that have more than 1 NA ---- 
All_batches <- All_batches %>% filter(NA_Count<2) %>% select(-NA_Count)

# ---- Save distance score list without AstroJ ---- 
All_batches  %>% arrange(Mean_DistanceScore_AllB_no_J) %>%
  write.csv(file = "SCR4_Confirmation_screening.csv")

# ---- Plot histogram of mean Distance Score for all batches ---- 
hist_DS_AllBatches <- All_batches %>% ggplot(aes(x=Mean_DistanceScore_AllB_no_J)) + 
  geom_histogram(fill= "#155F83E5") +
  geom_vline(xintercept = 0.7, linetype= "dashed", colour= "#8A9045E5", size=1)+
  ggtitle("All batches except AstroJ")+
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )
ggsave(file="Supp5B_SCR4_hist_DS_AllBatches_noJ.pdf", hist_DS_AllBatches, path = "../figs/", width = 4, height=3)

# ---- Create a new table with compounds appearing only once for those that are hits with different concentrations ---- 
Hit_List_individual_compounds <- All_batches %>% filter(Mean_DistanceScore_AllB_no_J <0.7) %>% ddply("Compound", nrow)

# ---- Add library of origin ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Compound_origin = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Rottlerin", "Tyrphostin A9", "JFD00244", "Methyl 2,5-dihydroxycinnamate", "PTACH")){
    Hit_List_individual_compounds$Compound_origin[i] <- "PKE library"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Carbendazim", "Ceritinib")){
    Hit_List_individual_compounds$Compound_origin[i] <- "Anti-cancer drug library"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Cefprozil (monohydrate)", "Chlorhexidine", "Clemastine (fumarate)", 
                                                      "Clomipramine (hydrochloride)", "Doxazosin (mesylate)", "Fenbendazole", "Flubendazole", 
                                                      "Flunarizine (dihydrochloride)", "Hexestrol", "Mebendazole", "Methylergometrine (maleate)", "Niclosamide", 
                                                      "Penicillin G benzathine (tetrahydrate)", "Perphenazine", "Salmeterol", "Ziprasidone (hydrochloride monohydrate)" )){
    Hit_List_individual_compounds$Compound_origin[i] <- "FDA-approved drug library"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("BMS-3", "BMS-5", "FRAX1036", "Ro 31-8220 (mesylate)", "SR-3306")){
    Hit_List_individual_compounds$Compound_origin[i] <- "Manually curated drugs"}
}

# ---- Add if tested in human ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Tested_in_human = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Carbendazim", "Cefprozil (monohydrate)", "Ceritinib", "Chlorhexidine", "Clemastine (fumarate)",
                                                      "Clomipramine (hydrochloride)", "Doxazosin (mesylate)", "Flubendazole", "Flunarizine (dihydrochloride)", "Hexestrol",
                                                      "Mebendazole", "Methylergometrine (maleate)", "Niclosamide", "Penicillin G benzathine (tetrahydrate)", "Perphenazine",
                                                      "Salmeterol", "Ziprasidone (hydrochloride monohydrate)")){
    Hit_List_individual_compounds$Tested_in_human[i] <- "Yes"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("BMS-3", "BMS-5", "Fenbendazole", "FRAX1036", "JFD00244", "Methyl 2,5-dihydroxycinnamate", "PTACH", "Ro 31-8220 (mesylate)",
                                                      "Rottlerin", "SR-3306", "Tyrphostin A9")){
    Hit_List_individual_compounds$Tested_in_human[i] <- "No"}
}

# ---- If tested in human, status of compounds ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Status = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Cefprozil (monohydrate)", "Ceritinib", "Chlorhexidine", "Clemastine (fumarate)", "Clomipramine (hydrochloride)", 
                                                      "Doxazosin (mesylate)", "Flubendazole", "Flunarizine (dihydrochloride)", "Mebendazole", "Methylergometrine (maleate)", "Niclosamide",
                                                      "Penicillin G benzathine (tetrahydrate)", "Perphenazine", "Salmeterol", "Ziprasidone (hydrochloride monohydrate)")){
    Hit_List_individual_compounds$Status[i] <- "Approved"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Carbendazim")){
    Hit_List_individual_compounds$Status[i] <- "Tested in clinical trials"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Hexestrol")){
    Hit_List_individual_compounds$Status[i] <- "Withdrawn"}
}

# ---- Add if expected to cross BBB ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Cross_BBB = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Carbendazim", "Ceritinib", "Clemastine (fumarate)", "Clomipramine (hydrochloride)", "Flunarizine (dihydrochloride)",
                                                      "Hexestrol", "Mebendazole", "Methylergometrine (maleate)", "Niclosamide", "Perphenazine", "Ziprasidone (hydrochloride monohydrate)")){
    Hit_List_individual_compounds$Cross_BBB[i] <- "Yes"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Cefprozil (monohydrate)", "Chlorhexidine", "Doxazosin (mesylate)", "Fenbendazole", "Flubendazole", "Penicillin G benzathine (tetrahydrate)",
                                                      "Salmeterol")){
    Hit_List_individual_compounds$Cross_BBB[i] <- "No"}
}

# ---- Add therapeutic class ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Therapeutic_class = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Flubendazole")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Metabolism"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Carbendazim", "Ceritinib")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Oncology"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Cefprozil (monohydrate)", "Chlorhexidine", "Fenbendazole", "Mebendazole", "Niclosamide", "Penicillin G benzathine (tetrahydrate)")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Infectiology"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Clemastine (fumarate)")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Allergology"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Clomipramine (hydrochloride)", "Flunarizine (dihydrochloride)", "Perphenazine", "Ziprasidone (hydrochloride monohydrate)")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Central Nervous System"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Doxazosin (mesylate)")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Cardiovascular"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Hexestrol")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Endocrinology"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Methylergometrine (maleate)")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Neuromuscular"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Salmeterol")){
    Hit_List_individual_compounds$Therapeutic_class[i] <- "Respiratory"}
}

# ---- Add family ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Family_of_drugs = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Fenbendazole", "Flubendazole", "Mebendazole", "Niclosamide")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Anti-parasitic"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Carbendazim")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Fungicide"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Cefprozil (monohydrate)", "Penicillin G benzathine (tetrahydrate)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antibiotic"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Ceritinib")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antineoplastic"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Chlorhexidine")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antiseptic"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Clemastine (fumarate)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antihistamine"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Clomipramine (hydrochloride)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antidepressant"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Doxazosin (mesylate)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antihypertensive"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Flunarizine (dihydrochloride)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Anti-migraine"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Hexestrol")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Synthetic estrogen"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Methylergometrine (maleate)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Alkaloid"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Perphenazine", "Ziprasidone (hydrochloride monohydrate)")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Antipsychotic"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Salmeterol")){
    Hit_List_individual_compounds$Family_of_drugs[i] <- "Anti-asthma"}
}

# ---- Add Targets ---- 
Hit_List_individual_compounds <- Hit_List_individual_compounds %>% mutate(Targets = NA)
for(i in 1:nrow(Hit_List_individual_compounds)){
  if(Hit_List_individual_compounds$Compound[i] %in% c("Albendazole", "Carbendazim", "Fenbendazole", "Flubendazole", "Mebendazole")){
    Hit_List_individual_compounds$Targets[i] <- "Micro-tubules"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("BMS-3", "BMS-5", "Ceritinib", "FRAX1036", "Ro 31-8220 (mesylate)", "SR-3306")){
    Hit_List_individual_compounds$Targets[i] <- "Kinase"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Cefprozil (monohydrate)", "Chlorhexidine", "Penicillin G benzathine (tetrahydrate)")){
    Hit_List_individual_compounds$Targets[i] <- "Bacteria cell wall"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Clomipramine (hydrochloride)", "Clemastine (fumarate)", "Doxazosin (mesylate)", "Methylergometrine (maleate)", "Perphenazine", "Salmeterol", "Ziprasidone (hydrochloride monohydrate)")){
    Hit_List_individual_compounds$Targets[i] <- "Neurotransmitter Receptor"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Flunarizine (dihydrochloride)")){
    Hit_List_individual_compounds$Targets[i] <- "Ion channel"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("JFD00244", "PTACH")){
  Hit_List_individual_compounds$Targets[i] <- "Epigenetic"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Methyl 2,5-dihydroxycinnamate", "Tyrphostin A9")){
    Hit_List_individual_compounds$Targets[i] <- "Cell surface receptor"}
  if(Hit_List_individual_compounds$Compound[i] %in% c("Niclosamide", "Hexestrol")){
    Hit_List_individual_compounds$Targets[i] <- "Transcription factor"}
}

# ---- Pie chart library of origins (all hits) ---- 
Hit_list_lib_origin <- Hit_List_individual_compounds %>% ddply("Compound_origin", nrow)
Hit_list_lib_origin <- Hit_list_lib_origin %>% mutate(Fraction_origin = (V1/nrow(Hit_List_individual_compounds))*100) %>%
  mutate(Fraction_origin = round(Fraction_origin, digits=1))
Hit_list_lib_origin <- Hit_list_lib_origin %>%
  arrange(desc(Compound_origin)) %>%
  mutate(lab.ypos=cumsum(Fraction_origin)-0.5*Fraction_origin)

Lib_origin_PieChart<- Hit_list_lib_origin %>% 
  ggplot(aes(x="", y=Fraction_origin, fill=Compound_origin))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal_ordered) +
  geom_text(aes(y = lab.ypos, label = V1), size=4, color = "black")+
  theme_void()+
  ggtitle("Library of origin")
ggsave(file="Fig2F_SCR4_PieChart_library_origin.pdf", Lib_origin_PieChart, path = "./results/", height = 4, width = 4)
rm(Hit_list_lib_origin, Lib_origin_PieChart)

# ---- Pie chart tested in human (all hits) ---- 
Hit_list_tested_human <- Hit_List_individual_compounds %>% ddply("Tested_in_human", nrow)
Hit_list_tested_human <- Hit_list_tested_human %>% mutate(Fraction_tested = (V1/nrow(Hit_List_individual_compounds))*100) %>%
  mutate(Fraction_tested = round(Fraction_tested, digits=1))
Hit_list_tested_human <- Hit_list_tested_human %>%
  arrange(desc(Tested_in_human)) %>%
  mutate(lab.ypos=cumsum(Fraction_tested)-0.5*Fraction_tested)
tested_human_PieChart<- Hit_list_tested_human %>% 
  ggplot(aes(x="", y=Fraction_tested, fill=Tested_in_human))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal_ordered_2) +
  geom_text(aes(y = lab.ypos, label = V1), size=4, color = "black")+
  theme_void()+
  ggtitle("Tested in human")
ggsave(file="Fig2F_SCR4_PieChart_tested_human.pdf", tested_human_PieChart, path = "./results/", height = 4, width = 4)
rm(Hit_list_tested_human, tested_human_PieChart)

# ---- Pie chart compounds target (all hits) ---- 
Hit_List_individual_compounds_targets <- Hit_List_individual_compounds %>% filter(!(is.na(Targets)))
Hit_list_target <- Hit_List_individual_compounds_targets %>% ddply("Targets", nrow)
Hit_list_target <- Hit_list_target %>% mutate(Fraction_target = (V1/nrow(Hit_List_individual_compounds_targets))*100) %>%
  mutate(Fraction_target = round(Fraction_target, digits=1))
Hit_list_target <- Hit_list_target %>%
  arrange(desc(Targets)) %>%
  mutate(lab.ypos=cumsum(Fraction_target)-0.5*Fraction_target)
Targets_PieChart <- Hit_list_target %>% 
  ggplot(aes(x="", y=Fraction_target, fill=Targets))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal) +
  geom_text(aes(y = lab.ypos, label = V1), size=4, color = "black")+
  theme_void()+
  ggtitle("Targets")
ggsave(file="Fig2F_SCR4_PieChart_Targets.pdf", Targets_PieChart, path = "./results/", height = 4, width = 4)
rm(Hit_List_individual_compounds_targets, Hit_list_target, Targets_PieChart)

# ---- Pie chart status of compounds (only hits tested in human) ---- 
Hit_List_individual_compounds_tested <- Hit_List_individual_compounds %>% filter(Tested_in_human == "Yes")
Hit_list_status <- Hit_List_individual_compounds_tested %>% ddply("Status", nrow)
Hit_list_status <- Hit_list_status %>% mutate(Fraction_status = (V1/nrow(Hit_List_individual_compounds_tested))*100) %>%
  mutate(Fraction_status = round(Fraction_status, digits=1))
Hit_list_status <- Hit_list_status %>%
  arrange(desc(Status)) %>%
  mutate(lab.ypos=cumsum(Fraction_status)-0.5*Fraction_status)
status_PieChart<- Hit_list_status %>% 
  ggplot(aes(x="", y=Fraction_status, fill=Status))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal_ordered_3) +
  geom_text(aes(y = lab.ypos, label = V1), size=4, color = "black")+
  theme_void()+
  ggtitle("Status of compound")
ggsave(file="Fig2F_SCR4_PieChart_status.pdf", status_PieChart, path = "./results/", height = 4, width = 4)
rm(Hit_list_status, status_PieChart)

# ---- Pie chart status of therapeutic class (only hits tested in human) ---- 
# Add therapeutic categories with those <2 classed as "others"
Hit_List_individual_compounds_tested <- Hit_List_individual_compounds_tested %>% mutate(Therapeutic_group = NA)
for(i in 1:nrow(Hit_List_individual_compounds_tested)){
  if(Hit_List_individual_compounds_tested$Therapeutic_class[i] %in% c("Allergology", "Cardiovascular", "Endocrinology", "Neuromuscular", "Respiratory")){
    Hit_List_individual_compounds_tested$Therapeutic_group[i] <- "Others"}
  if(Hit_List_individual_compounds_tested$Therapeutic_class[i] %in% c("Central Nervous System")){
    Hit_List_individual_compounds_tested$Therapeutic_group[i] <- "Central Nervous System"}
  if(Hit_List_individual_compounds_tested$Therapeutic_class[i] %in% c("Infectiology")){
    Hit_List_individual_compounds_tested$Therapeutic_group[i] <- "Infectiology"}
  if(Hit_List_individual_compounds_tested$Therapeutic_class[i] %in% c("Metabolism")){
    Hit_List_individual_compounds_tested$Therapeutic_group[i] <- "Metabolism"}
  if(Hit_List_individual_compounds_tested$Therapeutic_class[i] %in% c("Oncology")){
    Hit_List_individual_compounds_tested$Therapeutic_group[i] <- "Oncology"}
}
Hit_list_ther_group <- Hit_List_individual_compounds_tested %>% ddply("Therapeutic_group", nrow)
Hit_list_ther_group <- Hit_list_ther_group %>% mutate(Fraction_ther_group = (V1/nrow(Hit_List_individual_compounds_tested))*100) %>%
  mutate(Fraction_ther_group = round(Fraction_ther_group, digits=1))
Hit_list_ther_group <- Hit_list_ther_group %>%
  arrange(desc(Therapeutic_group)) %>%
  mutate(lab.ypos=cumsum(Fraction_ther_group)-0.5*Fraction_ther_group)
Ther_group_PieChart <- Hit_list_ther_group %>% 
  ggplot(aes(x="", y=Fraction_ther_group, fill=Therapeutic_group))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal_ordered_4) +
  geom_text(aes(y = lab.ypos, label = V1), size=4, color = "black")+
  theme_void()+
  ggtitle("Therapeutic class")
ggsave(file="Fig2F_SCR4_PieChart_therapeutic_group.pdf", Ther_group_PieChart, path = "./results/", height = 4, width = 4)
rm(Hit_list_ther_group, Ther_group_PieChart)

# ---- Pie chart expected to cross BBB (only hits tested in human) ---- 
Hit_list_BBB <- Hit_List_individual_compounds_tested %>% ddply("Cross_BBB", nrow)
Hit_list_BBB <- Hit_list_BBB %>% mutate(Fraction_BBB = (V1/nrow(Hit_List_individual_compounds_tested))*100) %>%
  mutate(Fraction_BBB = round(Fraction_BBB, digits=1))
Hit_list_BBB <- Hit_list_BBB %>%
  arrange(desc(Cross_BBB)) %>%
  mutate(lab.ypos=cumsum(Fraction_BBB)-0.5*Fraction_BBB)
BBB_PieChart<- Hit_list_BBB %>% 
  ggplot(aes(x="", y=Fraction_BBB, fill=Cross_BBB))+
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  scale_fill_manual(values = custom_pal_ordered_2) +
  geom_text(aes(y = lab.ypos, label = V1), size=4, color = "black")+
  theme_void()+
  ggtitle("Expected to cross BBB")
ggsave(file="Fig2F_SCR4_PieChart_BBB.pdf", BBB_PieChart, path = "./results/", height = 4, width = 4)
rm(Hit_list_status, status_PieChart)
