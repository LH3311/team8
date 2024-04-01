##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

# If you don't have BiocManager, here is the code to install it
# A lot of you probably already have this so you can skip
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none

# After installing all of its above dependencies, install ggpicrust2
install.packages("ggpicrust2")

#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(vegan)
library(dplyr)

# Set Working Directory
setwd("~/M - UBC/Courses/MICB 475/475_Group_Project/Aim_6")

#### Import files and preparing tables ####
#Importing the pathway PICrsut2
abundance_file <- "~/M - UBC/Courses/MICB 475/475_Group_Project/Aim_6/pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data <- as.data.frame(abundance_data)

# Get metadata with highlow
load("highlow_final.RData")
meta_hl <-  sample_data(pd_final)
metadata <-  data.frame(meta_hl)
metadata$ID = rownames(metadata)

#### Start with PD vs Control ####

#Remove NAs for your column of interest in this case subject        
metadata = metadata[!is.na(meta_hl$Disease),]


#Filtering the pathway table to only include samples that are in the filtered metadata        NOT NEEDED for Disease status, since there weren't any NA's in column
sample_names = metadata$'rownames'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data[, colSums(abundance_data != 0) > 0]

# Change the name of Column 1 in the Abundance Data
colnames(abundance_data_filtered)[1] <- 'pathway'

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadata = metadata[metadata$ID %in% abun_samples,] #making sure the filtered metadata only includes these samples

#Correct for individuals that are missing from metadata and not in the abundance matrix
metasamples = c("pathway",metadata$ID)
abundance_data_filtered = abundance_data_filtered[,colnames(abundance_data_filtered) %in% metasamples]


#### DESEq ####

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                        metadata = metadata, group = "Disease", daa_method = "DESeq2")

# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

feature_with_padjust_0.05 <- abundance_daa_results_df %>% filter(p_adjust < 0.05)

feature_with_padjust_0.01 <- abundance_daa_results_df %>% filter(p_adjust < 0.01)


## Using padjust < 0.05 ##

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_padjust_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_padjust_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_padjust_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
abundance_desc = abundance_desc[,-c(287:ncol(abundance_desc))] 

# Generate a heatmap
heatmap <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata, group = "Disease")
# save png file
ggsave(filename = "heatmap.png"
       , heatmap
       , height=40, width=40)

# Generate pathway PCA plot
PCA <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadata, group = "Disease")
# save png file
ggsave(filename = "PCA.png"
       , PCA
       , height=8, width=8)

## Using padjust < 0.05 ##

#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_padjust_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_padjust_0.05)

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_padjust_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description

#this line will change for each dataset. 34 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(285:ncol(abundance_desc))] 

# Generate a heatmap
heatmap2.1 <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata, group = "Disease")
# save png file
ggsave(filename = "heatmap2.1.png"
       , heatmap2.1
       , height=40, width=40)

# Generate pathway PCA plot
PCA2.1 <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadata, group = "Disease")
# save png file
ggsave(filename = "PCA2.1.png"
       , PCA
       , height=10, width=10)



# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "Disease")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(padj < 0.05)

# Filter by Log2fold change
sig_res_filtered = sig_res %>%
  filter(log2FoldChange <= -0.5 | log2FoldChange >= 0.5)

sig_res_filtered <- sig_res_filtered[order(sig_res_filtered$log2FoldChange),]
ggplot(data = sig_res_filtered, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")


#### for loop for all significant nutrients from Aim 1 ####

sig_nutrient <- c("Non_Starch_Polysaccharides_highlow", "Total_folate_highlow","Vitamin_C_highlow", "Vitamin_B2_highlow", "Fructose_highlow", "Coffe_drinker_highlow")

sig_pathways <- list()

for (nutrient in sig_nutrient) {
  
# Filter nutrient for PD

metadataPD = filter(metadata, !!sym(nutrient) == "pd_low"| !!sym(nutrient) == "pd_high")

  
# Filtering the pathway table to only include samples that are in the filtered metadata       
sample_names = metadataPD$'rownames'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering
  
#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data[, colSums(abundance_data != 0) > 0]
  
# Change the name of Column 1 in the Abundance Data
colnames(abundance_data_filtered)[1] <- 'pathway'
  
# Verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadataPD = metadataPD[metadataPD$ID %in% abun_samples,] #making sure the filtered metadata only includes these samples
  
  
# Correct for individuals that are missing from metadata and not in the abundance matrix
metasamples = c("pathway",metadataPD$ID)
abundance_data_filtered = abundance_data_filtered[,colnames(abundance_data_filtered) %in% metasamples]
  
# DESEq #
  
#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL
  
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                          metadata = metadataPD, group = nutrient, daa_method = "DESeq2")
  
# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                         daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE) 
  
# Filter p-values to only significant ones
sig_features <- abundance_daa_results_df %>% filter(p_adjust < 0.05)
  
#Changing the pathway column to description for the results 
feature_desc = inner_join(sig_features,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(sig_features)
  
#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% sig_features$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
abundance_desc = abundance_desc[,-c(287:ncol(abundance_desc))] 
  
# Generate a heatmap
heatmap <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadataPD, group = nutrient)
# save png file
ggsave(paste("/", nutrient, "/", "pd_heatmap.png", sep = ""), heatmap)
  
# Generate pathway PCA plot
pca <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadataPD, group = nutrient)
# save png file
ggsave(paste("/", nutrient, "/", "pd_pca.png", sep = ""), pca)
  
# Lead the function in
source("DESeq2_function.R")
  
# Run the function
res =  DEseq2_function(abundance_data_filtered, metadataPD, nutrient)
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)
  
# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(padj < 0.05)
  
sig_res <- sig_res[order(sig_res$log2FoldChange),]
log2foldchange -> ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
    geom_bar(stat = "identity")+ 
    theme_bw()+
    labs(x = "Log2FoldChange", y="Pathways")
  
ggsave(paste("/", nutrient, "/", "pd_log2foldchange.png", sep = ""), log2foldchange)

}  
  


# Filter nutrient for Control
  
metedataCON = filter(metadata, n %in% c("control_low", "control_high"))
  
# Filtering the pathway table to only include samples that are in the filtered metadata       
sample_names = metadataCON$'rownames'
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data[, colSums(abundance_data != 0) > 0]
  
# Change the name of Column 1 in the Abundance Data
colnames(abundance_data_filtered)[1] <- 'pathway'
  
# Verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
metadataCON = metadataCON[metadataCON$ID %in% abun_samples,] #making sure the filtered metadata only includes these samples

  
# Correct for individuals that are missing from metadata and not in the abundance matrix
metasamples = c("pathway",metadataCON$ID)
abundance_data_filtered = abundance_data_filtered[,colnames(abundance_data_filtered) %in% metasamples]
  
# DESEq #
  
#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL
  
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), 
                                          metadata = metadataCON, group = n, daa_method = "DESeq2")
  
# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                         daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE) 
  
# Filter p-values to only significant ones
sig_features <- abundance_daa_results_df %>% filter(p_adjust < 0.05)
  
#Changing the pathway column to description for the results 
feature_desc = inner_join(sig_features,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(sig_features)
  
#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% sig_features$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
abundance_desc = abundance_desc[,-c(287:ncol(abundance_desc))] 
  
# Generate a heatmap
heatmap <- pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadataCON, group = n)
# save png file
ggsave(paste("/", n, "/", "con_heatmap.png", sep = ""), heatmap)
  
# Generate pathway PCA plot
pca <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadataCON, group = n)
# save png file
ggsave(paste("/", n, "/", "con_pca.png", sep = ""), pca)
  
# Lead the function in
source("DESeq2_function.R")
  
# Run the function
res =  DEseq2_function(abundance_data_filtered, metadataCON, n)
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)
  
# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(padj < 0.05)
  
sig_res <- sig_res[order(sig_res$log2FoldChange),]
log2foldchange -> ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
    geom_bar(stat = "identity")+ 
    theme_bw()+
    labs(x = "Log2FoldChange", y="Pathways")
  
ggsave(paste("/", n, "/", "con_log2foldchange.png", sep = ""), log2foldchange)
  


