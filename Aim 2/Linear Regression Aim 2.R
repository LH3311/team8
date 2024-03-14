#Loading packages
#Code adapted from "Module15_TwoWayAnova_lm.R" by Dr. Evelyn Sun 
install.packages('car')
library(tidyverse)
library(phyloseq)
library(car)

#Loading rarified PD phyloseq object 
load("pd_rare.Rdata")

# Split into PD and controls 
pd_only_rare <- subset_samples(pd_rare, Disease == "PD")
control_rare <- subset_samples(pd_rare, Disease == "Control")


#Combining metadata and estimate_richness data

pd_only_wdiv <- data.frame(sample_data(pd_only_rare), estimate_richness(pd_only_rare))
control_wdiv <- data.frame(sample_data(control_rare), estimate_richness(control_rare))

#Defining predictor and response variable
pd_predictors <- colnames(pd_only_wdiv[, 38:98])
control_predictors <- colnames(control_wdiv[, 38:98])
response_variables <- c("Shannon", "Chao1", "ACE", "Simpson", "InvSimpson", "Fisher")

# Create formulas for linear model

pd_formulas <- lapply(response_variables, function(response_var) {
  as.formula(paste(response_var, "~", paste(pd_predictors, collapse = "+")))
})

control_formulas <- lapply(response_variables, function(response_var) {
  as.formula(paste(response_var, "~", paste(control_predictors, collapse = "+")))
})


# Fit linear models

pd_models <- lapply(pd_formulas, function(formula) {
  lm(formula, data = pd_only_wdiv)
})

control_models <- lapply(control_formulas, function(formula) {
  lm(formula, data = control_wdiv)
})

#For loop iterates through the two linear models to retrieve and print each alpha diversity metric 

for (i in seq_along(pd_models)) {
  cat("Summary of Model for Response Variable:", response_variables[i], "\n")
  print(summary(pd_models[[i]]))
  cat("\n")
}

for (i in seq_along(control_models)) {
  cat("Summary of Model for Response Variable:", response_variables[i], "\n")
  print(summary(control_models[[i]]))
  cat("\n")
}

# Creating empty list to store significant p values retrieved from for loops below
pd_significant_p_values <- list()
control_significant_p_values <- list()

# Access and extract significant p-values from each model
for (i in seq_along(pd_models)) {
  pd_model_summary <- summary(pd_models[[i]])
  pd_coefficients <- coef(pd_model_summary)
  pd_p_values <- pd_coefficients[, "Pr(>|t|)"]
  pd_significant_p_values[[response_variables[i]]] <- pd_p_values[pd_p_values < 0.05]
}

for (i in seq_along(control_models)) {
  control_model_summary <- summary(control_models[[i]])
  control_coefficients <- coef(control_model_summary)
  control_p_values <- control_coefficients[, "Pr(>|t|)"]
  control_significant_p_values[[response_variables[i]]] <- control_p_values[control_p_values < 0.05]
}

# Convert the list of significant p-values into a data frame
pd_significant_p_values_df <- do.call(rbind, pd_significant_p_values)
control_significant_p_values_df <- do.call(rbind, control_significant_p_values)

print(pd_significant_p_values_df)
print(control_significant_p_values_df)

#Testing av plots from car
for (i in seq_along(pd_models)) {
  cat("AV Plot for Response Variable:", response_variables[i], "\n")
  av_plot <- avPlots(pd_models[[i]], ask = FALSE)
  print(av_plot)
}

