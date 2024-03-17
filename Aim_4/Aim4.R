#install packages
install.packages("indicspecies") #install packages
library(tidyverse)
library(phyloseq)
library(indicspecies)

#load phyloseq object
load("../Downloads/highlow_final.RData")

#for loop of all nutrients
sig_nutrient <- c("Non_Starch_Polysaccharides_highlow", "Total_folate_highlow",
                  "Vitamin_C_highlow", "Vitamin_B2_highlow", "Fructose_highlow", "Coffe_drinker_highlow")

#group at genus level
highlow_genus <- tax_glom(pd_final, "Genus", NArm = FALSE)
#convert to relative abundance
highlow_genus_RA <- transform_sample_counts(highlow_genus, fun=function(x) x/sum(x))

#generate file for nutrient omitting NA
highlow_genus_RA_starchomit <- subset_samples(highlow_genus_RA, !is.na(sample_data(highlow_genus_RA)[["Non_Starch_Polysaccharides_highlow"]]))


isa_highlow_nonstarchpolysaccharides <- multipatt(t(otu_table(highlow_genus_RA_starchomit)), cluster = sample_data(highlow_genus_RA_starchomit)$`Non_Starch_Polysaccharides_highlow`)
summary(isa_highlow_nonstarchpolysaccharides)
taxtable <- tax_table(pd_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")
isa_highlow_nonstarchpolysaccharides$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()


