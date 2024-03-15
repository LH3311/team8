# Core Microbiome Analysis

#### Load Packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("highlow_final.RData")

#### Convert to relative abundance ####
pd_RA <- transform_sample_counts(pd_final, fun=function(x) x/sum(x))

#### for loop for all nutrient Venn Diagrams ####
sig_nutrient <- c("Non_Starch_Polysaccharides_highlow", "Total_folate_highlow",
                  "Vitamin_C_highlow", "Vitamin_B2_highlow", "Fructose_highlow", "Coffe_drinker_highlow")

core_member_tax_table_all <- list()

for (nutrient in sig_nutrient) {
  # Subset Data
  pd_stat_high <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "pd_high")
  pd_stat_low <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "pd_low")
  
  control_stat_high <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "control_high")
  control_stat_low <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "control_low")
  
  # Detection + Prevalence fractions
  d <- 0.01
  p <- 0.5
  
  pd_high_ASVs <- core_members(pd_stat_high, detection=d, prevalence=p)
  pd_low_ASVs <- core_members(pd_stat_low, detection=d, prevalence=p)
  
  control_high_ASVs <- core_members(control_stat_high, detection=d, prevalence=p)
  control_low_ASVs <- core_members(control_stat_low, detection=d, prevalence=p)
  
  # Venn diagram of all possible mixing of conditions
  con1 <- list(PD_High = pd_high_ASVs, PD_Low = pd_low_ASVs)
  con2 <- list(PD_High = pd_high_ASVs, Control_Low = control_low_ASVs)
  con3 <- list(Control_High = control_high_ASVs, Control_Low = control_low_ASVs)
  con4 <- list(Control_High = control_high_ASVs, PD_Low = pd_low_ASVs)
    
  venn1 = list(Venn(con1), "PD_highlow")
  venn2 = list(Venn(con2), "PD_high_Control_low")
  venn3 = list(Venn(con3), "Control_highlow")
  venn4 = list(Venn(con4), "Control_high_PD_low")
  
  venn_group <- list(venn1, venn2, venn3, venn4)
  
  for (v in venn_group) { 
    data = process_data(v[[1]], shape_id = "201")
    pd_venn <- plot_venn(data)
    
    ggsave(paste("Aim_3/", nutrient, "/", v[[2]] , ".png", sep = ""), pd_venn) # Printing Venn diagrams 
  }
  
  # Listing ASVs
  conditions <- list("pd_high_ASVs", "pd_low_ASVs", "control_high_ASVs", "control_low_ASVs")
  core_member_tax_tables <- list()
  for (cond in conditions) {
    tax_table_info <- prune_taxa(get(cond),pd_final) %>% tax_table()
    core_member_tax_tables[[cond]] <- list(key_nutrient = nutrient, tax_table = tax_table_info)
  }
  core_member_tax_table_all <- append(core_member_tax_table_all, core_member_tax_tables)
}

# Saving list of ASVs
save(core_member_tax_table_all, file="Aim_3/core_member_tax_table.RData")



#### Extra's to look at ASV relative abundance and membership ####

nutrient <- "Coffe_drinker_highlow" # enter nutrient of choice

prune_taxa(pd_high_ASVs,pd_RA) %>% # enter ASV group and subsets of interest
  plot_bar(fill="Order") +
  facet_wrap(.~get(nutrient), scales ="free")
  
con1 # PD_highlow
con2 # PD_high_Control_low
con3 # Control_highlow
con4 # Control_high_PD_low

