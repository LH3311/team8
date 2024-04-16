#import libraries 
library(tidyverse)
library(phyloseq)

# Load datasets
load("aim_1/pd_rare.RData")
load("aim_1/pd_final.RData")

# Split into PD and controls 
pd_only_rare <- subset_samples(pd_rare, Disease == "PD")
control_rare <- subset_samples(pd_rare, Disease == "Control")

# Split into PD and controls 
pd_only_final <- subset_samples(pd_final, Disease == "PD")
control_final <- subset_samples(pd_final, Disease == "Control")

sig_nutrient_names <- c("Non_Starch_Polysaccharides", "Total_folate",
                        "Vitamin_C", "Vitamin_B2", "Fructose")

############## High & Low Nutrients Non-Rarefied ##############


for (n in sig_nutrient_names) {
  # Choose the nutrient column and PD status column
  nutrient_column <- n      # Nutrient Column Name
  new_column_name <-  paste(nutrient_column,"_highlow", sep = "")       # New Column Name
  pd_status_column <- "Disease"                        # Disease Status
  
  # Calculate quartiles for the nutrient column
  nutrient_quartiles_pd <- quantile(sample_data(pd_only_final)[[nutrient_column]], c(0.25, 0.75))
  nutrient_quartiles_control <- quantile(sample_data(control_final)[[nutrient_column]], c(0.25, 0.75))
  
  # Create a new column based on PD status and nutrient quartiles
  sample_data(pd_final)[[new_column_name]] <- ifelse(
    sample_data(pd_final)[[nutrient_column]] <= nutrient_quartiles_pd[1] & sample_data(pd_final)$Disease == "PD",
    "pd_low",
    ifelse(
      sample_data(pd_final)[[nutrient_column]] <= nutrient_quartiles_control[1] & sample_data(pd_final)$Disease == "Control",
      "control_low",
      ifelse(
        sample_data(pd_final)[[nutrient_column]] >= nutrient_quartiles_pd[2] & sample_data(pd_final)$Disease == "PD",
        "pd_high",
        ifelse(
          sample_data(pd_final)[[nutrient_column]] >= nutrient_quartiles_control[2] & sample_data(pd_final)$Disease == "Control",
          "control_high",
          NA
        )
      )
    )
  )
}

###### Coffee Specifically

# Choose the nutrient column and PD status column
nutrient_column <- "Coffe_drinker"      # Nutrient Column Name
new_column_name <-  paste(nutrient_column,"_highlow", sep = "")       # New Column Name
pd_status_column <- "Disease"                        # Disease Status

# Calculate quartiles for the nutrient column
nutrient_quartiles_pd <- c("No", "Daily_multiple_cups")
nutrient_quartiles_control <- c("No", "Daily_multiple_cups")

# Create a new column based on PD status and nutrient quartiles
sample_data(pd_final)[[new_column_name]] <- ifelse(
  sample_data(pd_final)[[nutrient_column]] == nutrient_quartiles_pd[1] & sample_data(pd_final)$Disease == "PD",
  "pd_low",
  ifelse(
    sample_data(pd_final)[[nutrient_column]] == nutrient_quartiles_control[1] & sample_data(pd_final)$Disease == "Control",
    "control_low",
    ifelse(
      sample_data(pd_final)[[nutrient_column]] == nutrient_quartiles_pd[2] & sample_data(pd_final)$Disease == "PD",
      "pd_high",
      ifelse(
        sample_data(pd_final)[[nutrient_column]] == nutrient_quartiles_control[2] & sample_data(pd_final)$Disease == "Control",
        "control_high",
        NA
      )
    )
  )
)

# View metadata
meta_pd <-  sample_data(pd_final)
meta_pd <-  data.frame(meta_pd)

# Save file
save(pd_final, file="high_low/highlow_final.RData")

############## High & Low Nutrients Rarefied ##############

for (n in sig_nutrient_names) {
  # Choose the nutrient column and PD status column
  nutrient_column <- n      # Nutrient Column Name
  new_column_name <-  paste(nutrient_column,"_highlow", sep = "")       # New Column Name
  pd_status_column <- "Disease"                        # Disease Status
  
  # Calculate quartiles for the nutrient column
  nutrient_quartiles_pd <- quantile(sample_data(pd_only_rare)[[nutrient_column]], c(0.25, 0.75))
  nutrient_quartiles_control <- quantile(sample_data(control_rare)[[nutrient_column]], c(0.25, 0.75))
  
  # Create a new column based on PD status and nutrient quartiles
  sample_data(pd_rare)[[new_column_name]] <- ifelse(
    sample_data(pd_rare)[[nutrient_column]] <= nutrient_quartiles_pd[1] & sample_data(pd_rare)$Disease == "PD",
    "pd_low",
    ifelse(
      sample_data(pd_rare)[[nutrient_column]] <= nutrient_quartiles_control[1] & sample_data(pd_rare)$Disease == "Control",
      "control_low",
      ifelse(
        sample_data(pd_rare)[[nutrient_column]] >= nutrient_quartiles_pd[2] & sample_data(pd_rare)$Disease == "PD",
        "pd_high",
        ifelse(
          sample_data(pd_rare)[[nutrient_column]] >= nutrient_quartiles_control[2] & sample_data(pd_rare)$Disease == "Control",
          "control_high",
          NA
        )
      )
    )
  )
}

###### Coffee Specifically

# Choose the nutrient column and PD status column
nutrient_column <- "Coffe_drinker"      # Nutrient Column Name
new_column_name <-  paste(nutrient_column,"_highlow", sep = "")       # New Column Name
pd_status_column <- "Disease"                        # Disease Status

# Calculate quartiles for the nutrient column
nutrient_quartiles_pd <- c("No", "Daily_multiple_cups")
nutrient_quartiles_control <- c("No", "Daily_multiple_cups")

# Create a new column based on PD status and nutrient quartiles
sample_data(pd_rare)[[new_column_name]] <- ifelse(
  sample_data(pd_rare)[[nutrient_column]] == nutrient_quartiles_pd[1] & sample_data(pd_rare)$Disease == "PD",
  "pd_low",
  ifelse(
    sample_data(pd_rare)[[nutrient_column]] == nutrient_quartiles_control[1] & sample_data(pd_rare)$Disease == "Control",
    "control_low",
    ifelse(
      sample_data(pd_rare)[[nutrient_column]] == nutrient_quartiles_pd[2] & sample_data(pd_rare)$Disease == "PD",
      "pd_high",
      ifelse(
        sample_data(pd_rare)[[nutrient_column]] == nutrient_quartiles_control[2] & sample_data(pd_rare)$Disease == "Control",
        "control_high",
        NA
      )
    )
  )
)

# View metadata
meta_pd <-  sample_data(pd_rare)
meta_pd <-  data.frame(meta_pd)

# Save file
save(pd_rare, file="high_low/highlow_rare.RData")


