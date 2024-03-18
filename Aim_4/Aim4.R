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

#generate file for nutrient omitting NA for each nutrient
for (nutrient in sig_nutrient) { #for each nutrient in the sig_nutrient category
  subset_data <- subset_samples(highlow_genus_RA, !is.na(sample_data(highlow_genus_RA)[[nutrient]])) #remove all of the NA samples
  isa_subset <- multipatt(t(otu_table(subset_data)), cluster = sample_data(subset_data)[[nutrient]]) #running the isa analysis on the samples without NA
  taxtable <- tax_table(pd_final) %>% as.data.frame() %>% rownames_to_column(var="ASV") #extract row of taxa table from phyloseq object and convert to new row
  isa_subset_table <- isa_subset$sign %>% #save the table in a new object
    rownames_to_column(var="ASV") %>% #make the ASV rows column names in isa_subset_table
    left_join(taxtable) %>% #leftjoin and align
    filter(p.value<0.05) #filter from the asv table to remove non signficant
  write.csv(isa_subset_table, file = paste0(nutrient, "_filtered_table.csv"), row.names = FALSE) #save the new file as a csv
}


