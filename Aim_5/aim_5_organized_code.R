#### Load Packages ####
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(dplyr)

load("~/Documents/GitHub/team8/R/high_low/highlow_final.RData")

#convert to deseq, remove NAs
highlow_plus1 <- transform_sample_counts(pd_final, function(x) x+1)


sig_nutrients <- c("Non_Starch_Polysaccharides_highlow","Total_folate_highlow", 
                   "Vitamin_C_highlow", "Vitamin_B2_highlow", "Fructose_highlow", "Coffe_drinker_highlow")
  
### Settings, CHOOSE BEFORE RUNNING ANY OF THE FOLLOWING CODE ###
# Use for comparing high and low in PD
setting <- c("pd_high","pd_low", ": PD high vs PD low", "_PD_vol_plot.png")
# Use for comparing high and low in Control
setting <- c("control_high","control_low", ": Control high vs Control low", "_control_vol_plot.png")

#### Volcano plot ####

for (nutrient in sig_nutrients) {
  formula <- as.formula(paste("~", nutrient))
  highlow_deseq <- phyloseq_to_deseq2(prune_samples(!is.na(sample_data(highlow_plus1)[[nutrient]]) , highlow_plus1), formula)
  DESEQ <- DESeq(highlow_deseq)
  
  deseq_results <- results(DESEQ, tidy=TRUE, 
                        # right most value is your reference group
                        contrast = c(nutrient, setting[1],setting[2]))
  #View(deseq_results)
  #ggplot(deseq_results) +
  #  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
  
  ## Make variable to color by whether it is significant + large change
  title <- gsub("_", " ", nutrient)
  title <- gsub(" highlow", setting[3], title)
  
  PD_vol_plot <- deseq_results %>%
    mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
    ggplot() +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"))
    ggtitle(title)
  ## Save volcano plot
  ggsave(paste("Aim_5/", nutrient, setting[4], sep = ""),PD_vol_plot)
  
  print(paste("Completed", nutrient, setting[3], "Volcano plot"))
}

#### Bar Plot ####

for (nutrient in sig_nutrients) {
  formula <- as.formula(paste("~", nutrient))
  highlow_deseq <- phyloseq_to_deseq2(prune_samples(!is.na(sample_data(highlow_plus1)[[nutrient]]) , highlow_plus1), formula)
  DESEQ <- DESeq(highlow_deseq)
  
  deseq_results <- results(DESEQ, tidy=TRUE, 
                        # right most value is your reference group
                        contrast = c(nutrient, setting[1],setting[2]))
  
  ## To get table of results
  sigASVs <- deseq_results %>% 
    filter(padj<0.01 & abs(log2FoldChange)>2) %>%
    dplyr::rename(ASV=row)
  #View(sigASVs)
  ## Get only asv names
  sigASVs_vec <- sigASVs %>%
    pull(ASV)
  ## Prune phyloseq file
  title <- gsub("_", " ", nutrient)
  title <- gsub(" highlow", setting[3], title)
  
  highlow_deseq <- prune_taxa(sigASVs_vec,pd_final)
  sigASVs <- tax_table(highlow_deseq) %>% as.data.frame() %>%
    rownames_to_column(var="ASV") %>%
    right_join(sigASVs) %>%
    arrange(log2FoldChange) %>%
    mutate(Genus = make.unique(Genus)) %>%
    mutate(Genus = factor(Genus, levels=unique(Genus)))
  
  ## Generate barplot
  sigASVs_barplot <- ggplot(sigASVs) +
    geom_bar(aes(x=Genus, y=log2FoldChange, fill = Phylum), stat="identity")+
    geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    ggtitle(paste("Significant ASVs, ", title, sep = ""))
  sigASVs_barplot

  ## Save barplot
  ggsave(paste("Aim_5/sigASVs_", nutrient, setting[4], sep = ""), sigASVs_barplot)
  
  print(paste("Completed", nutrient, setting[3], "Bar plot"))
}

#### ASV list ####

# Complete list of ASVs
asv_list <- list()

for (nutrient in sig_nutrients) {
  formula <- as.formula(paste("~", nutrient))
  highlow_deseq <- phyloseq_to_deseq2(prune_samples(!is.na(sample_data(highlow_plus1)[[nutrient]]) , highlow_plus1), formula)
  DESEQ <- DESeq(highlow_deseq)
  
  deseq_results <- results(DESEQ, tidy=TRUE, 
                           # right most value is your reference group
                           contrast = c(nutrient, setting[1],setting[2]))
  
  ## To get table of results
  sigASVs <- deseq_results %>% 
    filter(padj<0.01 & abs(log2FoldChange)>2) %>%
    dplyr::rename(ASV=row)
  View(sigASVs)
  ## Get only asv names
  sigASVs_vec <- sigASVs %>%
    pull(ASV)
  ## Prune phyloseq file
  highlow_deseq <- prune_taxa(sigASVs_vec,pd_final)
  sigASVs <- tax_table(highlow_deseq) %>% as.data.frame() %>%
    rownames_to_column(var="ASV") %>%
    right_join(sigASVs) %>%
    arrange(log2FoldChange) %>%
    mutate(Genus = make.unique(Genus)) %>%
    mutate(Genus = factor(Genus, levels=unique(Genus)))
  
  ## Filter and Save ASV list
  top_5 <- sigASVs %>% 
    top_n(5, wt = log2FoldChange)
  lowest_5 <- sigASVs %>% 
    top_n(-5, wt = log2FoldChange)
  lowest_top_5 <- rbind(lowest_5, top_5)
  
  asv_nutrient <- list(key_nutrient = nutrient, condition = setting[3], ASVs = lowest_top_5)
  asv_list <- append(asv_list, asv_nutrient)
  
  print(paste("Completed", nutrient, setting[3], "ASV list append"))
}

# Saving list of ASVs
save(asv_list, file="Aim_5/asv_list.RData")

