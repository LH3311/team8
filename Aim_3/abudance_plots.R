# Raw Abundance

##########################################################

# Load Packages
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(microbiome)

load("highlow_final.RData")

# Specify the phylum, class, order, genus etc that you are interested in
taxa = "g__Agathobacter"
taxa_level = "Genus"

sig_taxa <- c("g__Bacteroides", "g__Agathobacter")

# Specify which nutrient to look in
nutrient <- "Fructose_highlow"
sig_nutrient <- c("Non_Starch_Polysaccharides_highlow", "Total_folate_highlow",
                  "Vitamin_C_highlow", "Vitamin_B2_highlow", "Fructose_highlow", "Coffe_drinker_highlow")

##########################################################

#Convert matrix to relative abundance
pd_RA <- transform_sample_counts(pd_final, fun=function(x) x/sum(x))
#Subset the phylsoeq to the taxa of interest
sub_beta_RA <- subset_taxa(pd_RA, Genus == taxa) ###################### MAKE SURE TAXA LEVEL MATCHES
#Collapse the ASVs of similar taxa rank (EX. Make all species of the same genus into one row)
sub_beta_order_RA <-tax_glom(sub_beta_RA, taxrank = taxa_level, NArm = FALSE)

#Extracting asv matrix and adding a sampleID column to merge with the metadata
asv_df = data.frame(otu_table(sub_beta_order_RA))
asv_df = data.frame(t(asv_df))
colnames(asv_df) = "relative_abundance"
asv_df$sampleID = rownames(asv_df)

#Extracting metadata and adding sampleID column to merge with asvs
metadata = data.frame(sample_data(sub_beta_order_RA))
metadata$sampleID = rownames(metadata)

#Merging ASVs and metadata, removing NA's from column of interest
asv_metadata_joined = inner_join(metadata, asv_df, by = "sampleID")
asv_metadata_joined = asv_metadata_joined[!is.na(asv_metadata_joined[[nutrient]]),]

#Plotting relative abundance as a barplot (violin plot also acceptable)
ggplot(asv_metadata_joined, aes(.data[[nutrient]], relative_abundance, fill = .data[[nutrient]]))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.title = element_text(size=rel(2), face = "bold"),
        axis.text = element_text(size=rel(1.4),face="bold", angle = 90),
        legend.position = "none")+
  labs(y="Relative Abundance (%)",x= "Treatments")


#Taking the average of relative abundance for all 4 groups
sum_df = asv_metadata_joined %>%
  group_by(Disease, .data[[nutrient]]) %>%
  dplyr::summarise(mean = mean(relative_abundance), sd = sd(relative_abundance))

ggplot(sum_df, aes(.data[[nutrient]], mean, fill = .data[[nutrient]]))+
  geom_errorbar(data = sum_df, aes(ymin = 0.005 ,ymax = mean+sd),width = 0.2, size = 1 )+
  geom_col()+
  theme_classic()+
  theme(axis.title = element_text(size=rel(2), face = "bold"),
        axis.text = element_text(size=rel(1.4),face="bold", angle = 90),
        legend.position = "none")+
  labs(y="Relative Abundance (%)",x= "Treatments")



#used for looking at taxa_table
#View(tax_table(pd_RA))

