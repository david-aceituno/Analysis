# MAIC comparing Ra223+ENZ vs Pluvicto
# Input: PEACE-3 IPD prepared in PEACE3_data_preparation_MAIC.R
# Output: MAIC-adjusted OS and PFS tables
# Author: David Aceituno
# Date: 2025-11-18

# 1. Packages used --------------------------------------------------------
pkgs <- c("tidyverse", "readxl", "here", "survival", "survminer", "haven")
lapply(pkgs, library, character.only = T)
source(here("R", "utils.R"))
km_directory <- here("data", "Pluvicto MAICs", "Reconstructed IPD")

# 2. Import data  --------------------------------------------------------
# PEACE IPD prepared in PEACE3_data_preparation_MAIC.R
peace3_ipd <- readRDS(here("data", "peace3_ipd_for_MAIC.RDS"))

# Import reconstructed IPD from VISION
comparator_ipd <- list()
comparator_ipd$OS <- read_csv(here(km_directory, "OS_VISION_PLU.csv"))
comparator_ipd$PFS <- read_csv(here(km_directory, "PFS_VISION_PLU.csv"))

# Load comparator baseline data
comparator_baseline <- as.data.frame(read_xlsx(here("data", "baseline_VISION.xlsx"))) #note: ECOG in VISION include ECOG 0-1
comparator_baseline <- comparator_baseline |> mutate(psa_log = log(psa))

# 3. Prepare data for MAIC --------------------------------------------------------
outcome_names <- c("OS", "PFS")
matching_variables <- c("age", "ecog", "gleason_8", "psa_log")

# Comparator name
comparator_names <- c("Pluvicto")

# MAIC model
km_fit <- cox_fit <- list()
maic_km_fit <- maic_cox_fit <- list()
maic_ipd <- list()

# Table to feed with with results
results_table <- matrix(nrow = 2, ncol = 7)
rownames(results_table) <- paste(rep(c("OS", "PFS"), length(comparator_names)), 
                                 rep(comparator_names, each = length(outcome_names)))
colnames(results_table) <- c("Naive median survival Ra223+ENZ", 
                             "Naive median survival Comparator",
                             "Naive HR",
                             "MAIC median survival Ra223+ENZ",
                             "MAIC median survival Comparator",
                             "MAIC HR",
                             "ESS")

comparator_name <- comparator_names[1]
comparator_study <- "VISION"

# 2. Run MAIC --------------------------------------------------------
for(outcome_name in outcome_names) {
  # Remove the placebo patients from peace3 IPD and combine with comparator data
  maic_ipd[[outcome_name]] <- peace3_ipd[[outcome_name]][peace3_ipd[[outcome_name]]$treatment == "Radium-223 + Enzalutamide", ]
  maic_ipd[[outcome_name]]  <- maic_ipd[[outcome_name]][, c("treatment", "time", "event", matching_variables)]
  
  # Harmonise names and change time to years
  #colnames(maic_ipd[[outcome_name]]) <- c("treatment", "time", "event", matching_variables)
  #maic_ipd[[outcome_name]]$time <- as.numeric(maic_ipd[[outcome_name]]$time) / 365.25
  
  # Combine with comparator data
  temp <- as.data.frame(comparator_ipd[[outcome_name]])
  temp <- temp[, c("time", "event")]
  temp$treatment <- comparator_name
  
  # Add empty data for the matching characteristics
  for(i in 1:length(matching_variables)) temp[[matching_variables[i]]] <- NA
  
  # Change time from months to years
  temp$time <- as.numeric(temp$time) / 12
  temp <- temp[, c("treatment", "time", "event", matching_variables)]
  maic_ipd[[outcome_name]] <- rbind(maic_ipd[[outcome_name]], temp)
  
  # Naive survival analyses
  km_fit[[outcome_name]] <- survfit(Surv(time, event) ~ treatment, data = maic_ipd[[outcome_name]])
  cox_fit[[outcome_name]] <- coxph(Surv(time, event) ~ treatment, data = maic_ipd[[outcome_name]])
  
  # KM plot
  
  # MAIC Reweighting
  # Number of Ra223+ENZ patients in peace3
  n_ipd_patients <- sum(maic_ipd[[outcome_name]]$treatment == "Radium-223 + Enzalutamide") 
  # Matching variables in peace3
  x_ipd <- maic_ipd[[outcome_name]][maic_ipd[[outcome_name]]$treatment == "Radium-223 + Enzalutamide",c(matching_variables)] 
  # Average matching variables in comparator
  x_comparator <- comparator_baseline[comparator_baseline$study == comparator_study & 
                                        comparator_baseline$trt == comparator_name,
                                      c(matching_variables)] 
  n_matching <- length(matching_variables)
  
  # Ensure it is a matrix so we can use matrix multiplication
  x_ipd <- as.matrix(x_ipd)
  
  # Set missing values in the IPD to their average
  for(k in 1:n_matching) {
    if(sum(is.na(x_ipd[, k])) > 0) x_ipd[is.na(x_ipd[, k]), k] <- mean(x_ipd[, k], na.rm = TRUE)
  }
  
  # center the IPD matching variables on comparator means but keep an uncentered version
  x_ipd_uncentered <- x_ipd
  
  # Requires all variables to be numeric
  for (k in 1:n_matching) {
    x_ipd[,k] <- x_ipd[,k] - x_comparator[,k]
  }
  
  # objective function to be minimized for weight estimation
  Q <- function(alpha , x_ipd) {
    return( sum( exp (x_ipd %*% alpha )))
  }
  
  alpha <- rep(1, n_matching) # arbitrary starting point for the optimiser
  
  # objective function minimized using BFGS
  Q_min <- optim(fn = Q, x_ipd = as.matrix(x_ipd), par=alpha, method ="BFGS")
  
  hat_alpha <- Q_min$par # finite solution is the logistic regression parameters
  log_hat_w <- rep(0, n_ipd_patients)
  for (k in 1:n_matching) {
    log_hat_w <- log_hat_w + hat_alpha[k] * x_ipd[,k]
  }
  
  # Estimated weights
  hat_w <- exp (log_hat_w) 
  
  # Estimated effective sample size
  ess <- sum (hat_w)^2 / sum (hat_w^2)
  
  baseline_table <- matrix(nrow = 3, ncol = n_matching)
  rownames(baseline_table) <- c(comparator_name, "Ra223+ENZ unweighted", "Ra223+ENZ weighted")
  colnames(baseline_table) <- matching_variables
  
  # Comparator means
  baseline_table[comparator_name, c(1:n_matching)] <- as.matrix(x_comparator)
  
  # Unweighted means
  baseline_table["Ra223+ENZ unweighted", ] <- colMeans(x_ipd_uncentered)
  
  # Weighted means
  baseline_table["Ra223+ENZ weighted", ] <- (hat_w %*% x_ipd_uncentered) / sum(hat_w)
  
  write.csv(baseline_table, here("results", paste0("baseline_table_", comparator_name, ".csv")))
  
  # Export a histogram of patient weights
  jpeg(here("results", paste0("weights_histogram_", comparator_name, "_", outcome_name, ".jpg")))
  hist(hat_w, breaks = 20, main = "Histogram of patient weights", xlab = "")
  dev.off()
  
  # Append weights = 1 for each comparator patient so they're unaffected
  maic_weights <- c(hat_w, rep(1, sum(maic_ipd[[outcome_name]]$treatment == comparator_name)))
  
  # MAIC weighted survival analysis
  maic_cox_fit[[outcome_name]] <- coxph(Surv(time, event) ~ treatment,
                                          robust = TRUE , weights = maic_weights, 
                                          data = maic_ipd[[outcome_name]] )
  
  maic_km_fit[[outcome_name]] <- survfit(Surv(time, event) ~ treatment, weights = maic_weights,
                                         data = maic_ipd[[outcome_name]])
  
  # Summarise the results with KM curves
  # KM unweighted (naive)
  jpeg(here("results", paste0("naive_km_", comparator_name, "_", outcome_name, ".jpg")))
  print(ggsurvplot(km_fit[[outcome_name]], data = maic_ipd[[outcome_name]], risk.table = TRUE, risk.table.height = 0.35, legend.labs =
                     c("Pluvicto", "Radium-223 + Enzalutamide"), xlab="Time in years", ylab = outcome_name,
                   linetype = c(1:2)))
  dev.off()
  
  # KM weighed by MAIC
  jpeg(here("results", paste0("maic_km_", comparator_name, "_", outcome_name, ".jpg")))
  print(ggsurvplot(maic_km_fit[[outcome_name]], data = maic_ipd[[outcome_name]], risk.table = TRUE, risk.table.height = 0.35, legend.labs =
                     c("Pluvicto", "Radium-223 + Enzalutamide"), xlab="Time in years", ylab = outcome_name,
                   linetype = c(1:2)))
  dev.off()
  
  # Summarise results in a table
  results_table[paste(outcome_name, comparator_name), "ESS"] <- format(ess, digits = 5)
  results_table[paste(outcome_name, comparator_name), "Naive median survival Ra223+ENZ"] <-
    format_continuous_result(12 * surv_median(km_fit[[outcome_name]])[2, c("median", "lower", "upper")])
  results_table[paste(outcome_name, comparator_name), "Naive median survival Comparator"] <-
    format_continuous_result(12 * surv_median(km_fit[[outcome_name]])[1, c("median", "lower", "upper")])
  
  results_table[paste(outcome_name, comparator_name), "Naive HR"] <-
    format_continuous_result(c(summary(cox_fit[[outcome_name]])$conf.int[1, c("exp(coef)", "lower .95", "upper .95")],
                               summary(cox_fit[[outcome_name]])$logtest["pvalue"]))
  
  results_table[paste(outcome_name, comparator_name), "MAIC median survival Ra223+ENZ"] <-
    format_continuous_result(12 * surv_median(maic_km_fit[[outcome_name]])[2, c("median", "lower", "upper")])
  results_table[paste(outcome_name, comparator_name), "MAIC median survival Comparator"] <-
    format_continuous_result(12 * surv_median(maic_km_fit[[outcome_name]])[1, c("median", "lower", "upper")])
  
  results_table[paste(outcome_name, comparator_name), "MAIC HR"] <-
    format_continuous_result(c(summary(maic_cox_fit[[outcome_name]])$conf.int[1, c("exp(coef)", "lower .95", "upper .95")],
                               summary(maic_cox_fit[[outcome_name]])$logtest["pvalue"]))
  
}

# 5. Export results --------------------------------------------------------
write.csv(results_table, here("results", "maic_results_VISION.csv"))  
