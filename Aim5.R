library(tidyverse)
library(phyloseq)
library(DESeq2)
library(dplyr)


load("~/Documents/GitHub/team8/R/high_low/highlow_final.RData")

###Transformation of Data
highlow_plus1 <- transform_sample_counts(pd_final, function(x) x+1)

#Remove all NAs from NSP highlow
highlow_NSPs <- prune_samples(!is.na(highlow_plus1@sam_data$Non_Starch_Polysaccharides_highlow), highlow_plus1)

#convert to deseq
highlow_deseq <- phyloseq_to_deseq2(highlow_NSPs, ~ Non_Starch_Polysaccharides_highlow)

####NSP DESeq
DESEQ_NSPs <- DESeq(highlow_deseq)
###NSP PD high vs low
#PD Volcano
NSPres_PD <- results(DESEQ_NSPs, tidy=TRUE, 
               #this will ensure that pd_low is your reference group
               contrast = c("Non_Starch_Polysaccharides_highlow","pd_high","pd_low"))
View(NSPres_PD)

ggplot(NSPres_PD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
NSPres_PD_vol_plot <- NSPres_PD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("Non Starch Polysaccharides: PD high vs PD low")
NSPres_PD_vol_plot

ggsave(filename="NSPPres_PD_vol_plot.png",NSPres_PD_vol_plot)

# To get table of results
sigASVs_NSP_PD <- NSPres_PD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_NSP_PD)
# Get only asv names
sigASVs_NSP_PD_vec <- sigASVs_NSP_PD %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq <- prune_taxa(sigASVs_NSP_PD_vec,pd_final)
sigASVs_NSP_PD <- tax_table(highlow_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_NSP_PD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_NSP_PD_barplot <- ggplot(sigASVs_NSP_PD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: NSPs PD high vs PD low")
sigASVs_NSP_PD_barplot

ggsave(filename="sigASVs_NSP_PD.png",sigASVs_NSP_PD_barplot)

#Control Volcano
NSPres_control <- results(DESEQ_NSPs, tidy=TRUE, 
                     #this will ensure that pd_low is your reference group
                     contrast = c("Non_Starch_Polysaccharides_highlow","control_high","control_low"))
View(NSPres_control)

ggplot(NSPres_control) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
NSPres_control_vol_plot <- NSPres_control %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Non Starch Polysaccharides: control high vs control low")
NSPres_control_vol_plot

ggsave(filename="NSPres_control_vol_plot.png",NSPres_control_vol_plot)

# To get table of results
sigASVs_NSP_control <- NSPres_control %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_NSP_control)
# Get only asv names
sigASVs_NSP_control_vec <- sigASVs_NSP_control %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq <- prune_taxa(sigASVs_NSP_control_vec,pd_final)
sigASVs_NSP_control <- tax_table(highlow_deseq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_NSP_control) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_NSP_control_barplot <- ggplot(sigASVs_NSP_control) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: NSPs control high vs control low")
sigASVs_NSP_control_barplot

ggsave(filename="sigASVs_NSP_control.png",sigASVs_NSP_control_barplot)

