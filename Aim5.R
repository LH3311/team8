library(tidyverse)
library(phyloseq)
library(DESeq2)
library(dplyr)


load("~/Documents/GitHub/team8/R/high_low/highlow_final.RData")

###Transformation of Data
highlow_plus1 <- transform_sample_counts(pd_final, function(x) x+1)

#convert to deseq, remove NAs
highlow_deseq_NSP <- phyloseq_to_deseq2(prune_samples(!is.na(highlow_plus1@sam_data$Non_Starch_Polysaccharides_highlow), highlow_plus1), ~ Non_Starch_Polysaccharides_highlow)


####NSP DESeq
DESEQ_NSPs <- DESeq(highlow_deseq_NSP)
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
highlow_deseq_NSP_PD <- prune_taxa(sigASVs_NSP_PD_vec,pd_final)
sigASVs_NSP_PD <- tax_table(highlow_deseq_NSP_PD) %>% as.data.frame() %>%
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
highlow_deseq_NSP_control <- prune_taxa(sigASVs_NSP_control_vec,pd_final)
sigASVs_NSP_control <- tax_table(highlow_deseq_NSP_control) %>% as.data.frame() %>%
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


#### Folate
#convert to deseq, remove NAs
highlow_deseq_TF <- phyloseq_to_deseq2(prune_samples(!is.na(highlow_plus1@sam_data$Total_folate_highlow), highlow_plus1), ~ Total_folate_highlow)

DESEQ_TF <- DESeq(highlow_deseq_TF)
###TF PD high vs low
#PD Volcano
TFres_PD <- results(DESEQ_TF, tidy=TRUE, 
                     #this will ensure that pd_low is your reference group
                     contrast = c("Total_folate_highlow","pd_high","pd_low"))
View(TFres_PD)

ggplot(TFres_PD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
TFres_PD_vol_plot <- TFres_PD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("Total Folate: PD high vs PD low")
TFres_PD_vol_plot

ggsave(filename="TFres_PD_vol_plot.png",TFres_PD_vol_plot)

# To get table of results
sigASVs_TF_PD <- TFres_PD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_TF_PD)
# Get only asv names
sigASVs_TF_PD_vec <- sigASVs_TF_PD %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_TF_PD <- prune_taxa(sigASVs_TF_PD_vec,pd_final)
sigASVs_TF_PD <- tax_table(highlow_deseq_TF_PD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_TF_PD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_TF_PD_barplot <- ggplot(sigASVs_TF_PD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: TF PD high vs PD low")
sigASVs_TF_PD_barplot

ggsave(filename="sigASVs_TF_PD.png",sigASVs_TF_PD_barplot)

#TF Control Volcano
TFres_control <- results(DESEQ_TF, tidy=TRUE, 
                          #this will ensure that pd_low is your reference group
                          contrast = c("Total_folate_highlow","control_high","control_low"))
View(TFres_control)

ggplot(TFres_control) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
TFres_control_vol_plot <- TFres_control %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Total Folate: control high vs control low")
TFres_control_vol_plot

ggsave(filename="TFres_control_vol_plot.png",TFres_control_vol_plot)

# To get table of results
sigASVs_TF_control <- TFres_control %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_TF_control)
# Get only asv names
sigASVs_TF_control_vec <- sigASVs_TF_control %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_TF_control <- prune_taxa(sigASVs_TF_control_vec,pd_final)
sigASVs_TF_control <- tax_table(highlow_deseq_TF_control) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_TF_control) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_TF_control_barplot <- ggplot(sigASVs_TF_control) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: TF control high vs control low")
sigASVs_TF_control_barplot

ggsave(filename="sigASVs_TF_control.png",sigASVs_TF_control_barplot)




#### Vitamin C
#convert to deseq, remove NAs
highlow_deseq_VitC <- phyloseq_to_deseq2(prune_samples(!is.na(highlow_plus1@sam_data$Vitamin_C_highlow), highlow_plus1), ~ Vitamin_C_highlow)

DESEQ_VitC <- DESeq(highlow_deseq_VitC)
###VitC PD high vs low
#PD Volcano
VitCres_PD <- results(DESEQ_VitC, tidy=TRUE, 
                    #this will ensure that pd_low is your reference group
                    contrast = c("Vitamin_C_highlow","pd_high","pd_low"))
View(VitCres_PD)

ggplot(VitCres_PD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
VitCres_PD_vol_plot <- VitCres_PD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("Vitamin C: PD high vs PD low")
VitCres_PD_vol_plot

ggsave(filename="VitCres_PD_vol_plot.png",VitCres_PD_vol_plot)

# To get table of results
sigASVs_VitC_PD <- VitCres_PD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_VitC_PD)
# Get only asv names
sigASVs_VitC_PD_vec <- sigASVs_VitC_PD %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_VitC_PD <- prune_taxa(sigASVs_VitC_PD_vec,pd_final)
sigASVs_VitC_PD <- tax_table(highlow_deseq_VitC_PD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_VitC_PD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_VitC_PD_barplot <- ggplot(sigASVs_VitC_PD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: VitC PD high vs PD low")
sigASVs_VitC_PD_barplot

ggsave(filename="sigASVs_VitC_PD.png",sigASVs_VitC_PD_barplot)

#VitC Control Volcano
VitCres_control <- results(DESEQ_VitC, tidy=TRUE, 
                         #this will ensure that pd_low is your reference group
                         contrast = c("Vitamin_C_highlow","control_high","control_low"))
View(VitCres_control)

ggplot(VitCres_control) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
VitCres_control_vol_plot <- VitCres_control %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Vitamin C: control high vs control low")
VitCres_control_vol_plot

ggsave(filename="VitCres_control_vol_plot.png",VitCres_control_vol_plot)

# To get table of results
sigASVs_VitC_control <- VitCres_control %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_VitC_control)
# Get only asv names
sigASVs_VitC_control_vec <- sigASVs_VitC_control %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_VitC_control <- prune_taxa(sigASVs_VitC_control_vec,pd_final)
sigASVs_VitC_control <- tax_table(highlow_deseq_VitC_control) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_VitC_control) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_VitC_control_barplot <- ggplot(sigASVs_VitC_control) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: VitC control high vs control low")
sigASVs_VitC_control_barplot

ggsave(filename="sigASVs_VitC_control.png",sigASVs_VitC_control_barplot)




#### Vitamin B2
#convert to deseq, remove NAs
highlow_deseq_VitB2 <- phyloseq_to_deseq2(prune_samples(!is.na(highlow_plus1@sam_data$Vitamin_B2_highlow), highlow_plus1), ~ Vitamin_B2_highlow)

DESEQ_VitB2 <- DESeq(highlow_deseq_VitB2)
###Vit B2 PD high vs low
#PD Volcano
VitB2res_PD <- results(DESEQ_VitB2, tidy=TRUE, 
                      #this will ensure that pd_low is your reference group
                      contrast = c("Vitamin_B2_highlow","pd_high","pd_low"))
View(VitB2res_PD)

ggplot(VitB2res_PD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
VitB2res_PD_vol_plot <- VitB2res_PD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("Vitamin B2: PD high vs PD low")
VitB2res_PD_vol_plot

ggsave(filename="VitB2res_PD_vol_plot.png",VitB2res_PD_vol_plot)

# To get table of results
sigASVs_VitB2_PD <- VitB2res_PD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_VitB2_PD)
# Get only asv names
sigASVs_VitB2_PD_vec <- sigASVs_VitB2_PD %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_VitB2_PD <- prune_taxa(sigASVs_VitB2_PD_vec,pd_final)
sigASVs_VitB2_PD <- tax_table(highlow_deseq_VitB2_PD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_VitB2_PD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_VitB2_PD_barplot <- ggplot(sigASVs_VitB2_PD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: VitB2 PD high vs PD low")
sigASVs_VitB2_PD_barplot

ggsave(filename="sigASVs_VitB2_PD.png",sigASVs_VitB2_PD_barplot)

#Control Volcano
VitB2res_control <- results(DESEQ_VitB2, tidy=TRUE, 
                           #this will ensure that pd_low is your reference group
                           contrast = c("Vitamin_B2_highlow","control_high","control_low"))
View(VitB2res_control)

ggplot(VitB2res_control) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
VitB2res_control_vol_plot <- VitB2res_control %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Vitamin B2: control high vs control low")
VitB2res_control_vol_plot

ggsave(filename="VitB2res_control_vol_plot.png",VitB2res_control_vol_plot)

# To get table of results
sigASVs_VitB2_control <- VitB2res_control %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_VitB2_control)
# Get only asv names
sigASVs_VitB2_control_vec <- sigASVs_VitB2_control %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_VitB2_control <- prune_taxa(sigASVs_VitB2_control_vec,pd_final)
sigASVs_VitB2_control <- tax_table(highlow_deseq_VitB2_control) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_VitB2_control) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_VitB2_control_barplot <- ggplot(sigASVs_VitB2_control) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: VitB2 control high vs control low")
sigASVs_VitB2_control_barplot

ggsave(filename="sigASVs_VitB2_control.png",sigASVs_VitB2_control_barplot)




#### Fructose
#convert to deseq, remove NAs
highlow_deseq_Fruc <- phyloseq_to_deseq2(prune_samples(!is.na(highlow_plus1@sam_data$Fructose_highlow), highlow_plus1), ~ Fructose_highlow)

DESEQ_Fruc <- DESeq(highlow_deseq_Fruc)
###Fruc PD high vs low
#PD Volcano
Frucres_PD <- results(DESEQ_Fruc, tidy=TRUE, 
                       #this will ensure that pd_low is your reference group
                       contrast = c("Fructose_highlow","pd_high","pd_low"))
View(Frucres_PD)

ggplot(Frucres_PD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
Frucres_PD_vol_plot <- Frucres_PD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("Fructose: PD high vs PD low")
Frucres_PD_vol_plot

ggsave(filename="Frucres_PD_vol_plot.png",Frucres_PD_vol_plot)

# To get table of results
sigASVs_Fruc_PD <- Frucres_PD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_Fruc_PD)
# Get only asv names
sigASVs_Fruc_PD_vec <- sigASVs_Fruc_PD %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_Fruc_PD <- prune_taxa(sigASVs_Fruc_PD_vec,pd_final)
sigASVs_Fruc_PD <- tax_table(highlow_deseq_Fruc_PD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_Fruc_PD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_Fruc_PD_barplot <- ggplot(sigASVs_Fruc_PD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: Fruc PD high vs PD low")
sigASVs_Fruc_PD_barplot

ggsave(filename="sigASVs_Fruc_PD.png",sigASVs_Fruc_PD_barplot)

#Control Volcano
Frucres_control <- results(DESEQ_Fruc, tidy=TRUE, 
                            #this will ensure that pd_low is your reference group
                            contrast = c("Fructose_highlow","control_high","control_low"))
View(Frucres_control)

ggplot(Frucres_control) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
Frucres_control_vol_plot <- Frucres_control %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Fructose: control high vs control low")
Frucres_control_vol_plot

ggsave(filename="Frucres_control_vol_plot.png",Frucres_control_vol_plot)

# To get table of results
sigASVs_Fruc_control <- Frucres_control %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_Fruc_control)
# Get only asv names
sigASVs_Fruc_control_vec <- sigASVs_Fruc_control %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_Fruc_control <- prune_taxa(sigASVs_Fruc_control_vec,pd_final)
sigASVs_Fruc_control <- tax_table(highlow_deseq_Fruc_control) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_Fruc_control) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_Fruc_control_barplot <- ggplot(sigASVs_Fruc_control) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: Fruc control high vs control low")
sigASVs_Fruc_control_barplot

ggsave(filename="sigASVs_Fruc_control.png",sigASVs_Fruc_control_barplot)




#### Coffee
#convert to deseq, remove NAs
highlow_deseq_COFE <- phyloseq_to_deseq2(prune_samples(!is.na(highlow_plus1@sam_data$Coffe_drinker_highlow), highlow_plus1), ~ Coffe_drinker_highlow)

DESEQ_COFE <- DESeq(highlow_deseq_COFE)
###COFE PD high vs low
#PD Volcano
COFEres_PD <- results(DESEQ_COFE, tidy=TRUE, 
                      #this will ensure that pd_low is your reference group
                      contrast = c("Coffe_drinker_highlow","pd_high","pd_low"))
View(COFEres_PD)

ggplot(COFEres_PD) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
COFEres_PD_vol_plot <- COFEres_PD %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  ggtitle("Coffee: PD high vs PD low")
COFEres_PD_vol_plot

ggsave(filename="COFEres_PD_vol_plot.png",COFEres_PD_vol_plot)

# To get table of results
sigASVs_COFE_PD <- COFEres_PD %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_COFE_PD)
# Get only asv names
sigASVs_COFE_PD_vec <- sigASVs_COFE_PD %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_COFE_PD <- prune_taxa(sigASVs_COFE_PD_vec,pd_final)
sigASVs_COFE_PD <- tax_table(highlow_deseq_COFE_PD) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_COFE_PD) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_COFE_PD_barplot <- ggplot(sigASVs_COFE_PD) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: COFE PD high vs PD low")
sigASVs_COFE_PD_barplot

ggsave(filename="sigASVs_COFE_PD.png",sigASVs_COFE_PD_barplot)

#Control Volcano
COFEres_control <- results(DESEQ_COFE, tidy=TRUE, 
                           #this will ensure that pd_low is your reference group
                           contrast = c("Coffe_drinker_highlow","control_high","control_low"))
View(COFEres_control)

ggplot(COFEres_control) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
COFEres_control_vol_plot <- COFEres_control %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  ggtitle("Coffee: control high vs control low")
COFEres_control_vol_plot

ggsave(filename="COFEres_control_vol_plot.png",COFEres_control_vol_plot)

# To get table of results
sigASVs_COFE_control <- COFEres_control %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_COFE_control)
# Get only asv names
sigASVs_COFE_control_vec <- sigASVs_COFE_control %>%
  pull(ASV)

# Prune phyloseq file
highlow_deseq_COFE_control <- prune_taxa(sigASVs_COFE_control_vec,pd_final)
sigASVs_COFE_control <- tax_table(highlow_deseq_COFE_control) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_COFE_control) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs_COFE_control_barplot <- ggplot(sigASVs_COFE_control) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  ggtitle("Significant ASVs: COFE control high vs control low")
sigASVs_COFE_control_barplot

ggsave(filename="sigASVs_COFE_control.png",sigASVs_COFE_control_barplot)

