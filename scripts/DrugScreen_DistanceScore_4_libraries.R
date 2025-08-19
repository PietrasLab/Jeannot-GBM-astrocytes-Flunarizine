# ================================================================
# Title: DrugScreen_DistanceScore_4_libraries
# Description: Code to produce the plots showing Distance scores for all primary drug screens and the confirmation screen
# (Fig 2B, 2E, Supp. Tables 3-5)
# ================================================================

# ---- Load libraries ----
library(tidyverse)
library(ggpubr)
library(ggsci)
library(wesanderson)

# ----  Import StratoMiner output files ---- 
SCR1_AstroD_Pears <- read.delim("./data/SCR1PKEAstroDAggData_Pearson_WithControls.txt")
SCR1_AstroF_Pears <- read.delim("./data/SCR1PKEAstroFAggData_Pearson_WithControls.txt")
SCR2_AstroD_Pears <- read.delim("./data/SCR2TargetMAstroDAggData_Pearson_WithControls.txt")
SCR2_AstroF_Pears <- read.delim("./data/SCR2TargetMAstroFAggData_Pearson_WithControls.txt")
SCR3_AstroD_Pears <- read.delim("./data/Scr3FDAAstroDAggData_Pearson_WithControls.txt")
SCR3_AstroF_Pears <- read.delim("./data/Scr3FDAAstroFAggData_Pearson_WithControls.txt")
SCR4_AstroG_Pears <- read.delim("./data/SCR4ConfScreenAstroG_withCTL.txt")
SCR4_AstroH_Pears <- read.delim("./data/SCR4ConfScreenAstroH_withCTL.txt")
SCR4_AstroI_Pears <- read.delim("./data/SCR4ConfScreenAstroI_withCTL.txt")
SCR4_AstroJ_Pears <- read.delim("./data/SCR4ConfScreenAstroJ_withCTL.txt")

# ---- Function to produce summary statistics in plots (mean and +/- sd) ---- 
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# ---- Prepare the color palettes ---- 
uchicago_colors <- pal_uchicago("default", alpha = 0.9)(5)
custom_order_colors <- uchicago_colors[c(4, 3, 5)] 
wes_colors <- wes_palette("IsleofDogs2", 5)
custom_order_colors_2 <- wes_colors[c(5, 1, 2, 3, 4)] 

# ---- Filter controls required for the analysis ---- 
SCR4_AstroG_Pears <- SCR4_AstroG_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 
SCR4_AstroH_Pears <- SCR4_AstroH_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 
SCR4_AstroI_Pears <- SCR4_AstroI_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 
SCR4_AstroJ_Pears <- SCR4_AstroJ_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 

# ----  Rename categories for all files  ---- 
SCR1_AstroD_Pears <- SCR1_AstroD_Pears %>%
  mutate(Category = case_when(Category == "SAMPLE" ~ "IR+ Compounds", Category == "NEG" ~ "IR+ DMSO", Category == "POS" ~ "IR- DMSO", TRUE ~ Category))
SCR1_AstroF_Pears <- SCR1_AstroF_Pears %>%
  mutate(Category = case_when(Category == "SAMPLE" ~ "IR+ Compounds", Category == "NEG" ~ "IR+ DMSO", Category == "POS" ~ "IR- DMSO", TRUE ~ Category))
SCR2_AstroD_Pears <- SCR2_AstroD_Pears %>%
  mutate(Category = case_when(Category == "SAMPLE" ~ "IR+ Compounds", Category == "NEG" ~ "IR+ DMSO", Category == "POS" ~ "IR- DMSO", TRUE ~ Category))
SCR2_AstroF_Pears <- SCR2_AstroF_Pears %>%
  mutate(Category = case_when(Category == "SAMPLE" ~ "IR+ Compounds", Category == "NEG" ~ "IR+ DMSO", Category == "POS" ~ "IR- DMSO", TRUE ~ Category))
SCR3_AstroD_Pears <- SCR3_AstroD_Pears %>%
  mutate(Category = case_when(Category == "SAMPLE" ~ "IR+ Compounds", Category == "NEG" ~ "IR+ DMSO", Category == "POS" ~ "IR- DMSO", TRUE ~ Category))
SCR3_AstroF_Pears <- SCR3_AstroF_Pears %>%
  mutate(Category = case_when(Category == "SAMPLE" ~ "IR+ Compounds", Category == "NEG" ~ "IR+ DMSO", Category == "POS" ~ "IR- DMSO", TRUE ~ Category))
SCR4_AstroG_Pears <- SCR4_AstroG_Pears %>%
  mutate(Category = case_when(Category == "Compounds" ~ "IR+ Compounds", Category == "NEG (DMSO IR)" ~ "IR+ DMSO", Category == "POS (DMSO)" ~ "IR- DMSO", TRUE ~ Category))
SCR4_AstroH_Pears <- SCR4_AstroH_Pears %>%
  mutate(Category = case_when(Category == "Compounds" ~ "IR+ Compounds", Category == "NEG (DMSO IR)" ~ "IR+ DMSO", Category == "POS (DMSO)" ~ "IR- DMSO", TRUE ~ Category))
SCR4_AstroI_Pears <- SCR4_AstroI_Pears %>%
  mutate(Category = case_when(Category == "Compounds" ~ "IR+ Compounds", Category == "NEG (DMSO IR)" ~ "IR+ DMSO", Category == "POS (DMSO)" ~ "IR- DMSO", TRUE ~ Category))
SCR4_AstroJ_Pears <- SCR4_AstroJ_Pears %>%
  mutate(Category = case_when(Category == "Compounds" ~ "IR+ Compounds", Category == "NEG (DMSO IR)" ~ "IR+ DMSO", Category == "POS (DMSO)" ~ "IR- DMSO", TRUE ~ Category))

#  ---- SCR1: Remove useless columns, add astrocyte batch and join tables  ---- 
SCR1_AstroD_Pears <- SCR1_AstroD_Pears  %>% mutate(Astrocyte = "AstroD") %>% 
  select(Compound, Category, Concentration, Astrocyte, DistanceScore, PCA01, PCA02, PCA03, PCA04, PCA05, PCA06, PCA07, PCA08, PCA09, PCA10)
SCR1_AstroF_Pears <- SCR1_AstroF_Pears %>% mutate(Astrocyte = "AstroF") %>%
  select(Compound, Category, Concentration, Astrocyte, DistanceScore, PCA01, PCA02, PCA03, PCA04, PCA05, PCA06, PCA07, PCA08, PCA09, PCA10)
SCR1 <- full_join(SCR1_AstroD_Pears, SCR1_AstroF_Pears)
SCR1 <- SCR1 %>% mutate(Library = "PKE library")
rm(SCR1_AstroD_Pears)
rm(SCR1_AstroF_Pears)

# ---- SCR2: Remove useless columns, add astrocyte batch and join tables ---- 
SCR2_AstroD_Pears <- SCR2_AstroD_Pears %>% mutate(Astrocyte = "AstroD") %>%
  select(Compound, Category, Concentration, Astrocyte, DistanceScore, PCA01, PCA02, PCA03, PCA04, PCA05, PCA06, PCA07)
SCR2_AstroF_Pears <- SCR2_AstroF_Pears %>% mutate(Astrocyte = "AstroF") %>%
  select(Compound, Category, Concentration, Astrocyte, DistanceScore, PCA01, PCA02, PCA03, PCA04, PCA05, PCA06, PCA07, PCA08, PCA09)
SCR2 <- full_join(SCR2_AstroD_Pears, SCR2_AstroF_Pears)
SCR2 <- SCR2 %>% mutate(Library = "Anti Cancer Drug Library")
rm(SCR2_AstroD_Pears)
rm(SCR2_AstroF_Pears)

# ---- SCR3: add information about Concentration, remove useless columns, add astrocyte batch and join tables ---- 
Concentration_D <- c(1:nrow(SCR3_AstroD_Pears))
Concentration_F <- c(1:nrow(SCR3_AstroF_Pears))
SCR3_AstroD_Pears <- SCR3_AstroD_Pears %>% mutate(Column = parse_number(Well))
SCR3_AstroF_Pears <- SCR3_AstroF_Pears %>% mutate(Column = parse_number(Well))
for(i in 1:nrow(SCR3_AstroD_Pears)){
  if(SCR3_AstroD_Pears$XPlate[i] %in% c("AstroD_IR_1uM_plate1", "AstroD_IR_1uM_plate2", "AstroD_IR_1uM_plate3", "AstroD_IR_1uM_plate4")){
    Concentration_D[i] <- "1uM"} 
  if(SCR3_AstroD_Pears$XPlate[i] %in% c("AstroD_IR_10uM_plate1", "AstroD_IR_10uM_plate2", "AstroD_IR_10uM_plate3", "AstroD_IR_10uM_plate4")){
    Concentration_D[i] <- "10uM"} 
  if(SCR3_AstroD_Pears$XPlate[i] == "AstroD_IR_1and10uM_plate5" & SCR3_AstroD_Pears$Column[i] %in% c(2:9)){
    Concentration_D[i] <- "1uM"} 
  if(SCR3_AstroD_Pears$XPlate[i] == "AstroD_IR_1and10uM_plate5" & SCR3_AstroD_Pears$Column[i] %in% c(12:19)){
    Concentration_D[i] <- "10uM"} 
}
for(i in 1:nrow(SCR3_AstroF_Pears)){
  if(SCR3_AstroF_Pears$XPlate[i] %in% c("AstroF_IR_1uM_plate1", "AstroF_IR_1uM_plate2", "AstroF_IR_1uM_plate3", "AstroF_IR_1uM_plate4")){
    Concentration_F[i] <- "1uM"} 
  if(SCR3_AstroF_Pears$XPlate[i] %in% c("AstroF_IR_10uM_plate1", "AstroF_IR_10uM_plate2", "AstroF_IR_10uM_plate3", "AstroF_IR_10uM_plate4")){
    Concentration_F[i] <- "10uM"} 
  if(SCR3_AstroF_Pears$XPlate[i] == "AstroF_IR_1and10uM_plate5" & SCR3_AstroF_Pears$Column[i] %in% c(2:9)){
    Concentration_F[i] <- "1uM"} 
  if(SCR3_AstroF_Pears$XPlate[i] == "AstroF_IR_1and10uM_plate5" & SCR3_AstroF_Pears$Column[i] %in% c(12:19)){
    Concentration_F[i] <- "10uM"} 
}
SCR3_AstroD_Pears <- SCR3_AstroD_Pears %>% mutate(Concentration = Concentration_D)
SCR3_AstroF_Pears <- SCR3_AstroF_Pears %>% mutate(Concentration = Concentration_F)
SCR3_AstroD_Pears <- SCR3_AstroD_Pears %>% mutate(Astrocyte = "AstroD") %>%
  select(Compound, Category, Concentration, Astrocyte, DistanceScore, PCA01, PCA02, PCA03, PCA04, PCA05, PCA06, PCA07, PCA08, PCA09, PCA10) 
SCR3_AstroF_Pears <- SCR3_AstroF_Pears %>% mutate(Astrocyte = "AstroF") %>%
  select(Compound, Category, Concentration, Astrocyte, DistanceScore, PCA01, PCA02, PCA03, PCA04, PCA05, PCA06, PCA07, PCA08, PCA09) 
SCR3 <- full_join(SCR3_AstroD_Pears, SCR3_AstroF_Pears)
SCR3 <- SCR3 %>% mutate(Library = "FDA Approved Drugs Library")
rm(SCR3_AstroD_Pears)
rm(SCR3_AstroF_Pears)

# ---- Save csv files for SCR1, SCR2, SCR3 ---- 
SCR1 %>% arrange(Compound) %>% write.csv(file = "SCR1_PKE_Library.csv")  
SCR2 %>% arrange(Compound) %>% write.csv(file = "SCR2_AntiCancer_Drug_Library.csv")  
SCR3 %>% arrange(Compound) %>% write.csv(file = "SCR3_FDA_Approved_Compounds_Library.csv")  

# ---- Join all libraries except SCR4 ---- 
pre_All_Libraries <- full_join(SCR1, SCR2)
All_librairies <- full_join(pre_All_Libraries, SCR3)
rm(pre_All_Libraries)

# ---- Remove useless columns, add astrocyte batch and join tables ---- 
SCR4_AstroG_Pears <- SCR4_AstroG_Pears %>% select(DistanceScore, Compound, Category, Concentration) %>% 
  mutate(Astrocyte = "AstroG")
SCR4_AstroH_Pears <- SCR4_AstroH_Pears %>% select(DistanceScore, Compound, Category, Concentration) %>% 
  mutate(Astrocyte = "AstroH")
SCR4_AstroI_Pears <- SCR4_AstroI_Pears %>% select(DistanceScore, Compound, Category, Concentration) %>% 
  mutate(Astrocyte = "AstroI")
SCR4_AstroJ_Pears <- SCR4_AstroJ_Pears %>% select(DistanceScore, Compound, Category, Concentration) %>% 
  mutate(Astrocyte = "AstroJ")
pre_SCR4_1 <- full_join(SCR4_AstroG_Pears, SCR4_AstroH_Pears)
pre_SCR4_2 <- full_join(pre_SCR4_1, SCR4_AstroI_Pears)
SCR4 <- full_join(pre_SCR4_2, SCR4_AstroJ_Pears)
rm(SCR4_AstroG_Pears)
rm(SCR4_AstroH_Pears)
rm(SCR4_AstroI_Pears)
rm(SCR4_AstroJ_Pears)
rm(pre_SCR4_1)
rm(pre_SCR4_2)

#  ---- Make a BoxPlot with faceting Astrocytes/Libraries ---- 
Boxplot_AllLibraries  <- All_librairies %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  mutate(Library = as.factor(Library)) %>%
  mutate(Library = fct_relevel(Library, "PKE library", "Anti Cancer Drug Library", "FDA Approved Drugs Library")) %>%
  ggplot(aes(Category, DistanceScore)) +
  geom_jitter(aes(color=Category), position = position_jitter(0.2), size=0.5) +
  stat_summary(fun.data=data_summary, col = "grey40") +
  facet_grid(Astrocyte~Library) +
  ylab("\n Distance Score (AU)") + #le \n permet d'éloigner le titre du graph#
  theme_bw() + 
  geom_hline(yintercept = 0.7, linetype= "dashed", colour= "grey30", size=1) +
  geom_hline(yintercept = 1, linetype= "dotted", size=0.8, colour= "grey30") +
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10, angle=45, hjust = 1),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  scale_color_manual(values = custom_order_colors)
ggsave(file="Fig2B_DistanceScore_All_Libraries.pdf", Boxplot_AllLibraries, path = "./results/", bg = "transparent", height = 4, width = 6)

# ---- Save the .csv files with distance scores for each screen ----
All_librairies  %>% filter(Library == "PKE library") %>%
  write.csv(file = "./results/SuppTable3_SCR1_PKE_Library.csv")
All_librairies  %>% filter(Library == "Anti Cancer Drug Library") %>%
  write.csv(file = "./results/SuppTable4_SCR2_AntiCancer_Drug_Library.csv")
All_librairies  %>% filter(Library == "FDA Approved Drugs Library") %>%
  write.csv(file = "./results/SuppTable5_SCR3_FDA_Approved_Compounds_Library.csv")

# Plot Distance score confirmation screening SCR4
BoxPlot_SCR4  <- SCR4 %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  mutate(Concentration = fct_relevel(Concentration, "", "0.25µM", "1µM", "3µM", "10µM")) %>%
  ggplot(aes(Category, DistanceScore)) +
  geom_jitter(aes(color=Concentration), position = position_jitter(0.2), size=1) +
  stat_summary(fun.data=data_summary, col = "grey40") +
  ggtitle("Confirmation screening") +
  facet_grid(.~Astrocyte) +
  ylab("\n Distance Score (AU)") + #le \n permet d'éloigner le titre du graph#
  theme_bw() +    
  geom_hline(yintercept = 0.7, linetype= "dashed", colour= "grey30", linewidth=1) +
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    axis.text.x = element_text(angle=45, hjust = 1),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  scale_color_manual(values = custom_order_colors_2)

ggsave(file="Fig2E_DistanceScore_Confirmation_Screening.pdf", BoxPlot_SCR4, path = "./results/", bg = "transparent", height = 3, width = 8)








