# Rob Cloke
# March 2, 2024


#Adapted from Christopher Lee's script_for_multivariable_analysis

# This script is designed to rank multiple variables in terms of their contribution to microbial diversity using a PERMANOVA and looking at the
# R- squared value and pvalue. The R-squared statistic will indicate how much any one variable contributes to driving diversity.

#set your working directory
setwd("C:/Users/cloke/Documents/M - UBC/Courses/MICB 475/475_Group_Project/Aim_1")


#import libraries 
library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(cowplot)

# Load rarafied dataset
load("pd_rare.RData")

# Split into PD and controls 
pd_only_rare <- subset_samples(pd_rare, Disease == "PD")
control_rare <- subset_samples(pd_rare, Disease == "Control")



#Take out the ASV matrix from the phyloseq object and convert it to a dataframe.
ASV_table_all <-  otu_table(pd_rare)
ASV_table_all <-  data.frame(ASV_table_all)
ASV_table_allt <- t(ASV_table_all)               # transpose so Sample ID's are the rownames

ASV_table_pd <-  otu_table(pd_only_rare)
ASV_table_pd <-  data.frame(ASV_table_pd)
ASV_table_pdt <- t(ASV_table_pd)               # transpose so Sample ID's are the rownames

ASV_table_control <-  otu_table(control_rare)
ASV_table_control <-  data.frame(ASV_table_control)
ASV_table_controlt <- t(ASV_table_control)               # transpose so Sample ID's are the rownames
  
#Take out the metadata from the phyloseq object and convert it to a dataframe.
meta_all <-  sample_data(pd_rare)
meta_all <-  data.frame(meta_all)

meta_pd <-  sample_data(pd_only_rare)
meta_pd <-  data.frame(meta_pd)

meta_control <-  sample_data(control_rare)
meta_control <-  data.frame(meta_control)

#### Entire Population ####

#Filter the metadata to remove columns that are not related to nutrients
#
meta_all_nutrients <- meta_all[,c(38:98)]


#Without filtering the metadata, the following line captures all column headers.
nutrients_a = colnames(meta_all_nutrients)

#Build an empty list that will be filled up by the loop
adonis.res_a = list()   

#Create a loop to go over each variable in the metadata.
for (i in 1:61){                                #COMPLETE the for loop argument. You need to to loop through as many variables that are present in "nutrients". Use a range, IE (1:?)
  print(i)                                                 #Printing i to keep track of loop progress
  print(nutrients_a[[i]])
  meta_df_a <- data.frame(meta_all_nutrients[!is.na(meta_all_nutrients[, i]), ])    #Remove the rows in metadata that contain missing data for the i'th variable
  
  samples = rownames(meta_df_a)                              #Create a vector of all the samples in the metadata after removing NA's
  ASV_mat <- ASV_table_allt[match(rownames(meta_df_a), rownames(ASV_table_allt)),]                                              #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
                                #This is important because we need the number of individuals represented in the ASV table and metadata to be the same.
  
  # I did this earlier ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
  dis = vegdist(ASV_mat, method = "bray",diag =T, upper = T) #Create the distance matrix based on a bray-curtis dissimilarity model
  adonis.res_a[[i]] = vegan::adonis2(as.formula(paste("dis~",nutrients_a[[i]],sep = "")), data = meta_df_a) #This line runs the PERMANOVA test and puts it into the empty list we made above.
}

#Create an empty table to import the R-squared and pvalue.
result_a = matrix(NA, nrow = length(nutrients_a), ncol =2)

#Loop though each variable and generate a table that can be plotted.
for (i in 1:61){ #This for loop argument will be the same as on line 60
  
  result_a[i,1] = adonis.res_a[[i]][1,3] #Grab the R-squared statistic
  result_a[i,2] = adonis.res_a[[i]][1,5]#Grab the pvalue

}

rownames(result_a) <-  c(nutrients_a)                           #Convert the rowmanes to variables 
colnames(result_a) <-  c("R2", "Pvalue")                      #Change the column names TO "R2" AND "Pvalue"
result_a <-  data.frame(result_a, stringsAsFactors = F)           #Convert it to a data.frame (easiest to work with when plotting)
result_a$Padjust <-  p.adjust(result_a$Pvalue, method = "fdr")    #Generate an adjusted pvalue to correct for the probability of false positives
result_a$Nutrient <-  c(nutrients_a) 
View(result_a)

save(result_a, file="result_a.RData")

###############################PLOTTING

#filter the results table to only include significant variables with a pvalue<0.05

result_filtered_a_p <- result_a %>% filter(Pvalue < 0.05)
result_filtered_a_pad <- result_a %>% filter(Padjust < 0.05)

# Visualize Nutrients with Significant P-values < 0.05

results_a_pvalue <- ggplot(data = result_filtered_a_p, aes(x = reorder(Nutrient, R2, decreasing = FALSE),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() + 
  labs(y = "Adonis R2", x = "Nutrient", title = "Nutrients with P-Value < 0.05 in Entire Population")
results_a_pvalue

# save png file
ggsave(filename = "results_a_pvalue.png"
       , results_a_pvalue
       , height=4, width=4)

# Visualize Nutrients with Significant adjusted P-values < 0.05

results_a_padjust <- ggplot(data = result_filtered_a_pad, aes(x = reorder(Nutrient, R2, decreasing = FALSE),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() + 
  labs(y = "Adonis R2", x = "Nutrient", title = "Nutrients with Adjusted P < 0.05 in Entire Population")
results_a_padjust

# save png
ggsave(filename = "results_a_pvadjust.png"
       , results_a_padjust
       , height=4, width=4)

#### PD Patrients ####

#Filter the metadata to remove columns that are not related to nutrients
#
meta_pd_nutrients <- meta_pd[,c(38:98)]


#Without filtering the metadata, the following line captures all column headers.
nutrients_pd = colnames(meta_pd_nutrients)

adonis.res_pd = list()   #Build an empty list that will be filled up by the loop

#Create a loop to go over each variable in the metadata.
for (i in 1:61){                                #COMPLETE the for loop argument. You need to to loop through as many variables that are present in "nutrients". Use a range, IE (1:?)
  print(i)                                                 #Printing i to keep track of loop progress
  print(nutrients_pd[[i]])
  meta_df_pd <- data.frame(meta_pd_nutrients[!is.na(meta_pd_nutrients[, i]), ])    #Remove the rows in metadata that contain missing data for the i'th variable
  
  samples = rownames(meta_df_pd)                              #Create a vector of all the samples in the metadata after removing NA's
  ASV_matp <- ASV_table_pdt[match(rownames(meta_df_pd), rownames(ASV_table_pdt)),]                                              #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
  #This is important because we need the number of individuals represented in the ASV table and metadata to be the same.
  
  # I did this earlier ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
  dis = vegdist(ASV_matp, method = "bray",diag =T, upper = T) #Create the distance matrix based on a bray-curtis dissimilarity model
  adonis.res_pd[[i]] = vegan::adonis2(as.formula(paste("dis~",nutrients_pd[[i]],sep = "")), data = meta_df_pd) #This line runs the PERMANOVA test and puts it into the empty list we made above.
}

#Create an empty table to import the R-squared and pvalue.
result_pd = matrix(NA, nrow = length(nutrients_pd), ncol =2)

#Loop though each variable and generate a table that can be plotted.
for (i in 1:61){ #This for loop argument will be the same as on line 60
  
  result_pd[i,1] = adonis.res_pd[[i]][1,3] #Grab the R-squared statistic
  result_pd[i,2] = adonis.res_pd[[i]][1,5]#Grab the pvalue
  
}

rownames(result_pd) <-  c(nutrients_pd)                           #Convert the rowmanes to variables 
colnames(result_pd) <-  c("R2", "Pvalue")                      #Change the column names TO "R2" AND "Pvalue"
result_pd <-  data.frame(result_pd, stringsAsFactors = F)           #Convert it to a data.frame (easiest to work with when plotting)
result_pd$Padjust <-  p.adjust(result_pd$Pvalue, method = "fdr")    #Generate an adjusted pvalue to correct for the probability of false positives
result_pd$Nutrient <-  c(nutrients_pd) 
View(result_pd)

save(result_pd, file="result_pd.RData")

###############################PLOTTING

#Try and generate the plot yourself using ggplot.
#Here is a skeleton of the plotting code. I wrote "FILLOUT" where information needs to be added.

#Also, filter the results table to only include significant variables with a pvalue<0.05

result_filtered_pd_p <- result_pd %>% filter(Pvalue < 0.05)
result_filtered_pd_pad <- result_pd %>% filter(Padjust < 0.05)

# Visualize Nutrients with Significant P-values < 0.05
results_pd_pvalue <- ggplot(data = result_filtered_pd_p, aes(x = reorder(Nutrient, R2, decreasing = FALSE),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  labs(y = "Adonis R2", x = "Nutrient", title = "Nutrients with P-Value < 0.05 in PD Patients")
results_pd_pvalue

# save png file
ggsave(filename = "results_pd_pvalue.png"
       , results_pd_pvalue
       , height=4, width=4)

# Visualize Nutrients with Significant P-values < 0.05
results_pd_padjust <- ggplot(data = result_filtered_pd_pad, aes(x = reorder(Nutrient,  R2),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() + 
  labs(y = "Adonis R2", x = "Nutrient", title = "Nutrients with Adjusted P < 0.05 in PD Patients")
results_pd_padjust 

# save png file
ggsave(filename = "results_pd_padjust.png"
       , results_pd_padjust
       , height=4, width=4)


#### Controls ####

#Filter the metadata to remove columns that are not related to nutrients

meta_c_nutrients <- meta_control[,c(38:98)]


#Without filtering the metadata, the following line captures all column headers.
nutrients_c = colnames(meta_c_nutrients)

adonis.res_c = list()   #Build an empty list that will be filled up by the loop

#Create a loop to go over each variable in the metadata.
for (i in 1:61){                                #COMPLETE the for loop argument. You need to to loop through as many variables that are present in "nutrients". Use a range, IE (1:?)
  print(i)                                                 #Printing i to keep track of loop progress
  print(nutrients_c[[i]])
  meta_df_c <- data.frame(meta_c_nutrients[!is.na(meta_c_nutrients[, i]), ])    #Remove the rows in metadata that contain missing data for the i'th variable
  
  samples = rownames(meta_df_c)                              #Create a vector of all the samples in the metadata after removing NA's
  ASV_matc <- ASV_table_controlt[match(rownames(meta_df_c), rownames(ASV_table_controlt)),]                                              #Filter the ASV_table to remove the same individuals that were filtered out for NA values 
  #This is important because we need the number of individuals represented in the ASV table and metadata to be the same.
  
  # I did this earlier ASV_mat = t(ASV_mat +1) # Transpose the ASV matrix so that the sample IDs become rownames. This is required for  adonis2()
  dis = vegdist(ASV_matc, method = "bray",diag =T, upper = T) #Create the distance matrix based on a bray-curtis dissimilarity model
  adonis.res_c[[i]] = vegan::adonis2(as.formula(paste("dis~",nutrients_c[[i]],sep = "")), data = meta_df_c) #This line runs the PERMANOVA test and puts it into the empty list we made above.
}

#Create an empty table to import the R-squared and pvalue.
result_c = matrix(NA, nrow = length(nutrients_c), ncol =2)

#Loop though each variable and generate a table that can be plotted.
for (i in 1:61){ #This for loop argument will be the same as on line 60
  
  result_c[i,1] = adonis.res_c[[i]][1,3] #Grab the R-squared statistic
  result_c[i,2] = adonis.res_c[[i]][1,5]#Grab the pvalue
  
}

rownames(result_c) <-  c(nutrients_c)                           #Convert the rowmanes to variables 
colnames(result_c) <-  c("R2", "Pvalue")                      #Change the column names TO "R2" AND "Pvalue"
result_c <-  data.frame(result_c, stringsAsFactors = F)           #Convert it to a data.frame (easiest to work with when plotting)
result_c$Padjust <-  p.adjust(result_c$Pvalue, method = "fdr")    #Generate an adjusted pvalue to correct for the probability of false positives
result_c$Nutrient <-  c(nutrients_c) 
View(result_c)

save(result_c, file="result_c.RData")

###############################PLOTTING

# filter the results table to only include significant variables with a pvalue<0.05
result_filtered_c_p <- result_c %>% filter(Pvalue < 0.05)
result_filtered_c_pad <- result_c %>% filter(Padjust < 0.05)

# Visualize Nutrients with Significant P-values < 0.05
results_c_pvalue <- ggplot(data = result_filtered_c_p, aes(x = reorder(Nutrient, R2, decreasing = FALSE),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() +
  labs(y = "Adonis R2", x = "Nutrient", title = "Nutrients with P-Value < 0.05 in Control Particiapants")
results_c_pvalue

# save png file
ggsave(filename = "results_c_pvalue.png"
       , results_c_pvalue
       , height=4, width=4)

# Visualize Nutrients with Significant adjusted P-values < 0.05
results_c_padjust <- ggplot(data = result_filtered_c_pad, aes(x = reorder(Nutrient, R2, decreasing = FALSE),y=R2)) +
  geom_bar(stat='identity') +
  coord_flip() + 
  labs(y = "Adonis R2", x = "Nutrient", title = "Nutrients with Adjusted P < 0.05 in Control Participants")
results_c_padjust 

# save png file
ggsave(filename = "results_c_padjust.png"
       , results_c_padjust
       , height=4, width=4)

# Plot Entire Population, PD Patients and Controls in side-by-side graphs
p_value_comparison <- plot_grid(results_a_pvalue, results_pd_pvalue, results_c_pvalue, labels = "AUTO")

ggsave(filename = "p_value_comparison.png"
       , p_value_comparison 
       , height=10, width=10)

padjust_comparison <- plot_grid(results_a_padjust, results_pd_padjust, results_c_padjust, labels = "AUTO")

ggsave(filename = "padjust_comparison.png"
       , padjust_comparison 
       , height=10, width=10)
