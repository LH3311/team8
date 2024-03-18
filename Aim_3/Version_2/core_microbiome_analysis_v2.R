# Core Microbiome Analysis

#### Load Packages ####
library(tidyverse)
library(ggplot2)
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

  # list with the core_members() results for all 4 groups.
  nutrient_x <- list(PD_High = pd_high_ASVs, PD_Low = pd_low_ASVs, Control_High = control_high_ASVs, Control_Low = control_low_ASVs )

  # Venn diagram plotting
  pd_venn <- ggVennDiagram(nutrient_x,
                label = "count",
                edge_size = 0.5,
                category.names = c("PD High","PD Low","Control High","Control Low"))+
    ggplot2::scale_fill_gradient(low="white",high = "darkorchid4")+
    scale_x_continuous(expand = expansion(mult = .2))
  
  ggsave(paste("Aim_3/", nutrient, "_venn_diagram.png", sep = ""), pd_venn) # Printing Venn diagrams 
}

#### for loop for getting unique ASVs from all groups  ####

complete_asv_table <- data.frame()

for (nutrient in sig_nutrient) {
  # Subset Data
  pd_stat_high <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "pd_high")
  pd_stat_low <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "pd_low")
  
  control_stat_high <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "control_high")
  control_stat_low <- subset_samples(pd_RA, sample_data(pd_RA)[[nutrient]] == "control_low")
  
  # Detection + Prevalence fractions
  d <- 0.01
  p <- 0.25
  
  pd_high_ASVs <- core_members(pd_stat_high, detection=d, prevalence=p)
  pd_low_ASVs <- core_members(pd_stat_low, detection=d, prevalence=p)
  
  control_high_ASVs <- core_members(control_stat_high, detection=d, prevalence=p)
  control_low_ASVs <- core_members(control_stat_low, detection=d, prevalence=p)
  
  # Getting unique ASVs from all groups 
  total_res = append(pd_high_ASVs,pd_low_ASVs)
  total_res = append(total_res,control_high_ASVs)
  total_res = append(total_res,control_low_ASVs)
  total_res = data.frame(total_res)
  names = as.matrix(unique(total_res))

  # Making empty lists
  shared_all_asv = list()
  pd_high_asv = list()
  pd_low_asv = list()
  
  # Extracting ASVs of interest
  for (i in names){
    if (i %in% pd_high_ASVs & !i %in% pd_low_ASVs & !i %in% control_high_ASVs & !i %in% control_low_ASVs) { #PD high only
      pd_high_asv[[i]] = i
    } else if (!i %in% pd_high_ASVs & i %in% pd_low_ASVs & !i %in% control_high_ASVs & !i %in% control_low_ASVs) { # PD low only
      pd_low_asv[[i]] = i
    } else if (i %in% pd_high_ASVs & i %in% pd_low_ASVs & i %in% control_high_ASVs & i %in% control_low_ASVs) { #Shared in all groups
      shared_all_asv[[i]] = i
    }
  }
  all_asv <- list(
    shared_all_asv = shared_all_asv,
    pd_high_asv = pd_high_asv,
    pd_low_asv = pd_low_asv
  )
  
  # Turning ASVs of interest into table with taxa info
  taxa_info = data.frame(tax_table(pd_RA))
  taxa_info$ASV = rownames(taxa_info)
  asv_table = stack(all_asv)
  colnames(asv_table) = c("ASV", "Exclusive ASVs")
  asv_table$Nutrient = nutrient
  asv_table$dp = paste(d,p)
  rownames(asv_table) = NULL
  asv_table$ASV = as.character(asv_table$ASV)
  shared_taxa = inner_join(taxa_info, asv_table, by = "ASV")
  
  #Adding into complete table
  complete_asv_table <- rbind(complete_asv_table, shared_taxa)
}

write.csv(complete_asv_table, "Aim_3/complete_asv_table.csv")
save(complete_asv_table, file="Aim_3/complete_asv_table.RData")
  