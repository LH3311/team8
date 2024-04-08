#Loading packages
#Code adapted from "Module15_TwoWayAnova_lm.R" by Dr. Evelyn Sun 
install.packages('car')
install.packages('broom')
library(tidyverse)
library(phyloseq)
library(broom)
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


#Above code modified to instead store all the calculated p_values and their respective response variable
pd_p_values_table <- data.frame(Response_Variable = character(),
                                P_Value = numeric(),
                                stringsAsFactors = FALSE)

control_p_values_table <- data.frame(Response_Variable = character(),
                                P_Value = numeric(),
                                stringsAsFactors = FALSE)

# Loop through each model and extract p-values
for (i in seq_along(pd_models)) {
  model_summary <- summary(pd_models[[i]])
  pd_p_values <- coef(model_summary)[, "Pr(>|t|)"]  
  
  response_variable <- response_variables[i]
  pd_p_values_table <- rbind(pd_p_values_table, data.frame(Response_Variable = response_variable,
                                                           P_Value = pd_p_values))
  
  cat("Summary of Model for Response Variable:", response_variable, "\n")
  print(model_summary)
  cat("\n")
}

for (i in seq_along(control_models)) {
  model_summary <- summary(control_models[[i]])
  control_p_values <- coef(model_summary)[, "Pr(>|t|)"]  
  
  response_variable <- response_variables[i]
  control_p_values_table <- rbind(control_p_values_table, data.frame(Response_Variable = response_variable,
                                                           P_Value = control_p_values))
  
  cat("Summary of Model for Response Variable:", response_variable, "\n")
  print(model_summary)
  cat("\n")
}



# Creating empty list to store significant p values retrieved from for loops below
pd_significant_p_values <- data.frame()
control_significant_p_values <- data.frame()

# Access and extract significant p-values from each model
for (i in seq_along(pd_models)) {
  pd_model_summary <- summary(pd_models[[i]])
  pd_coefficients <- as.data.frame(coef(pd_model_summary))
  pd_coefficients$predictors = rownames(pd_coefficients)
  
  pd_pval_filtered = pd_coefficients %>%
    filter(pd_coefficients$`Pr(>|t|)` < 0.05)
  
  if (dim(pd_pval_filtered)[1] != 0 ){
    pd_pval_filtered$response = response_variables[[i]]
    pd_significant_p_values <- rbind(pd_significant_p_values,pd_pval_filtered)
  } else{}
}

for (i in seq_along(control_models)) {
  #i=4
  control_model_summary <- summary(control_models[[i]])
  control_coefficients <- as.data.frame(coef(control_model_summary))
  control_coefficients$predictors = rownames(control_coefficients)
  
  control_pval_filtered = control_coefficients %>%
    filter(control_coefficients$`Pr(>|t|)` < 0.05)
  
  if (dim(control_pval_filtered)[1] != 0 ){
    control_pval_filtered$response = response_variables[[i]]
    control_significant_p_values <- rbind(control_significant_p_values,control_pval_filtered)
  } else{}
  
  
}


#Plotting Significant Non_alcoholic_bevs against ACE, Fisher, Chao1 

ace_gg <- ggplot(control_wdiv, aes(x = Non_alcoholic_bevs, y = ACE)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  xlab('Non-Alcoholic Beverages (g)') +
  ylab('ACE')

chao1_gg <- ggplot(control_wdiv, aes(x = Non_alcoholic_bevs, y = Chao1)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") + 
  xlab('Non-Alcoholic Beverages (g)') +
  ylab('Chao1') +
  theme(axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 12))

fisher_gg <- ggplot(control_wdiv, aes(x = Non_alcoholic_bevs, y = Fisher)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  xlab('Non-Alcoholic Beverages (g)') +
  ylab('Fisher')

#Plotting Significant Beta_carotene against Simpson
simpson_gg <- ggplot(control_wdiv, aes(x = Beta_carotene, y = Simpson)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  xlab('Beta Carotene (μg)') + 
  ylab('Simpson') +
  theme(axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 12))

ggsave("Non Alcoholic Beverages Against ACE final.png", plot = ace_gg, width = 5, height = 5)
ggsave("Non Alcoholic Beverages Against Chao1 final2.png", plot = chao1_gg, width = 5, height = 5)
ggsave("Non Alcoholic Beverages Against Fisher final.png", plot = fisher_gg, width = 5, height =5)
ggsave("Beta Carotene Against Simpson final2.png", plot = simpson_gg, width = 5, height = 5)


#Testing av plots from car
for (i in seq_along(pd_models)) {
  cat("AV Plot for Response Variable:", response_variables[i], "\n")
  av_plot <- avPlots(pd_models[[i]], ask = FALSE)
  print(av_plot)
}



