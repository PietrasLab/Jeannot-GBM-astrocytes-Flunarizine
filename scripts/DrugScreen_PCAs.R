# ================================================================
# Title: DrugScreen_DistanceScore_4_libraries
# Description: Code to produce the PCA plots showing for all primary drug screens and the confirmation screen
# (Fig Supp. 4B-D, Fig. Supp. 5A)
# ================================================================

# ---- Load libraries ----
library(tidyverse)
library(ggsci)

# ---- Import StratoMiner output files ---- 
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

# ---- Prepare the color palette ----
uchicago_colors <- pal_uchicago("default", alpha = 0.7)(5)
custom_order_colors <- uchicago_colors[c(4, 3, 5)] 

# ---- Keep only useful controls for SCR4 ----
SCR4_AstroG_Pears <- SCR4_AstroG_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 
SCR4_AstroH_Pears <- SCR4_AstroH_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 
SCR4_AstroI_Pears <- SCR4_AstroI_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 
SCR4_AstroJ_Pears <- SCR4_AstroJ_Pears %>% filter(Category %in% c("Compounds", "NEG (DMSO IR)", "POS (DMSO)")) 

# ---- Rename categories for all files ----
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

# ---- Plot the PCAs using StratoMiner values ----
SCR1_AstroD_PCA <- SCR1_AstroD_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("PKE library - AstroD") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig4B_PCA_SCR1_AstroD.pdf", SCR1_AstroD_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR1_AstroF_PCA <- SCR1_AstroF_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("PKE library - AstroF") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig4B_PCA_SCR1_AstroF.pdf", SCR1_AstroF_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR2_AstroD_PCA <- SCR2_AstroD_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("Anti cancer drugs library - AstroD") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig4C_PCA_SCR2_AstroD.pdf", SCR2_AstroD_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR2_AstroF_PCA <- SCR2_AstroF_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("Anti cancer drugs library - AstroF") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig4C_PCA_SCR2_AstroF.pdf", SCR2_AstroF_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR3_AstroD_PCA <- SCR3_AstroD_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category), size=0.4) +
  ggtitle("FDA approved compounds library - AstroD") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig4D_PCA_SCR3_AstroD.pdf", SCR3_AstroD_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR3_AstroF_PCA <- SCR3_AstroF_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category), size=0.4) +
  ggtitle("FDA approved compounds library  - AstroF") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig4D_PCA_SCR3_AstroF.pdf", SCR3_AstroF_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR4_AstroG_PCA <- SCR4_AstroG_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("Confirmation screening - AstroG") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig5A_PCA_SCR4_AstroG.pdf", SCR4_AstroG_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR4_AstroH_PCA <- SCR4_AstroH_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("Confirmation screening - AstroH") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig5A_PCA_SCR4_AstroH.pdf", SCR4_AstroH_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR4_AstroI_PCA <- SCR4_AstroI_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("Confirmation screening - AstroI") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig5A_PCA_SCR4_AstroI.pdf", SCR4_AstroI_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)

SCR4_AstroJ_PCA <- SCR4_AstroJ_Pears %>% mutate(Category = fct_relevel(Category, "IR- DMSO", "IR+ DMSO", "IR+ Compounds")) %>%
  ggplot(aes(PCA01, PCA02)) + 
  geom_point(aes(col=Category)) +
  ggtitle("Confirmation screening - AstroJ") +
  theme_bw() + 
  theme(
    plot.title = element_text(size=12),
    axis.title = element_text(size=10, face="bold"),
    axis.text = element_text(size=10),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )+
  scale_color_manual(values = custom_order_colors)
ggsave(file="SuppFig5A_PCA_SCR4_AstroJ.pdf", SCR4_AstroJ_PCA, path = "./results/", bg = "transparent", height = 3, width = 5)



