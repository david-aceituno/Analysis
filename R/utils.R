# Utils functions

# Function to extract fit statistics from a BUGS object
extract_fit <- function(bugs_fit){
  result_fit <- c(bugs_fit$DIC,
                  bugs_fit$summary["totresdev",][c(1,3,7)])
  return(result_fit)
}

# Function to format results in CI style
format_results <- function(x, n_digits = 3) {
  return(paste0(format(x[1], digits = n_digits, nsmall = n_digits),
                " (", format(x[2], digits = n_digits, nsmall = n_digits),
                ", ", format(x[3], digits = n_digits, nsmall = n_digits), ")"))
}

# Function to extract summary statistics using CI style

summary.stats<-function(x,n.digits=2,med=FALSE)
{
  if(med){
    return(paste0(round(median(x),digits=n.digits)," (",round(quantile(x,probs=c(0.025)),digits=n.digits),", ",round(quantile(x,probs=c(0.975)),digits=n.digits),")"))
  }else{
    return(paste0(round(mean(x),digits=n.digits)," (",round(quantile(x,probs=c(0.025)),digits=n.digits),", ",round(quantile(x,probs=c(0.975)),digits=n.digits),")"))
  }
}

# Generate all pairs of powers of fractional polynomials
generate_fp_so_powers <- function(fp_powers) {
  fp_so_powers <- as.data.frame(matrix(nrow = sum(1:length(fp_powers)),
                                       ncol = 2))
  colnames(fp_so_powers) <- c("P1", "P2")
  fp_so_powers[, "P1"] <- rep(fp_powers, times = seq(from = length(fp_powers), to = 1, by = -1))
  for(i in 1:length(fp_powers)) {
    fp_so_powers[fp_so_powers[, "P1"] == fp_powers[i], "P2"] <- 
      fp_powers[i:length(fp_powers)]
  }
  return(fp_so_powers)
}

# Time transformations for second order FP
fp_transform_times <- function(times, powers = c(-1, 0.5)) 
{
  nobs <- length(times)
  npoly <- length(powers)
  timetrans <- matrix(0, nrow = nobs, ncol = npoly)
  timetrans1 <- ifelse(powers[1] != rep(0, nobs), times^powers[1], log(times))
  timetrans[, 1] <- timetrans1
  if (npoly >= 2) {
    for (i in 2:npoly) {
      if (powers[i] == powers[(i - 1)])
        timetrans2 <- log(times) * timetrans1
      else timetrans2 <- ifelse(powers[i] != rep(0, nobs), times^powers[i],
                                log(times))
      timetrans[, i] <- timetrans2
      timetrans1 <- timetrans2
    }
  }
  return(timetrans)
}

# Function to generate initial values for FP models
fp_inits <- function(t, aux){
  mh <- muhaz::muhaz(t, max.time = max(t))
  lhdat <- data.frame(loghaz = log(mh$haz.est),
                      time = mh$est.grid)
  lhdat <- na.omit(lhdat)
  lhdat <- lhdat[is.finite(lhdat$loghaz),]
  X <- flexsurv:::bfp(lhdat$time, powers=aux$powers)
  coef(lm(loghaz ~ X, data=lhdat))
}

# Format a vector of mean, 2.5% CI and 97.5% CI
format_results <- function(x, ndigits = 2) {
  paste0(format(x[1], digits = ndigits, nsmall = ndigits), " (",
         format(x[2], digits = ndigits, nsmall = ndigits), ", ",
         format(x[3], digits = ndigits, nsmall = ndigits), ")")
}

# Function to summarise relative effects from multinma
# Designed for survival outcomes
format_relative_effects <- function(multinma_rel, hazard_scale = TRUE,
                                    invert = TRUE,
                                    treatment_names = NULL) {
  multinma_rel <- as.data.frame(multinma_rel$summary)
  formatted_results <- matrix(NA, nrow = nrow(multinma_rel))
  for(i in 1:nrow(multinma_rel)) {
    if(invert){
      if(hazard_scale) {
        # Exponentiate the log hazard ratios
        formatted_results[i, ] <- 
          format_results(exp(-multinma_rel[i, c("mean", "97.5%", "2.5%")])  )
      } else {
        formatted_results[i, ] <- 
          format_results(-multinma_rel[i, c("mean", "97.5%", "2.5%")])
      }
    } else {
      if(hazard_scale) {
        # Exponentiate the log hazard ratios
        formatted_results[i, ] <- 
          format_results(exp( multinma_rel[i, c("mean", "2.5%", "97.5%")])  )
      } else {
        formatted_results[i, ] <- 
          format_results(multinma_rel[i, c("mean", "2.5%", "97.5%")])
      }
    }
  }
  colnames(formatted_results) <- c("Mean (95% CrI)")
  # Use treatment names if provided.
  if(!is.null(treatment_names)) {
    rownames(formatted_results) <- treatment_names
  } else {
    rownames(formatted_results) <- multinma_rel[, "parameter"]
  }
  return(formatted_results)
}

extract_surv_fit <- function(surv_object){
  summary_surv <- cbind(surv_object$res.t, surv_object$knots, surv_object$cov)
  summary_surv <- as_tibble(summary_surv)
  names(summary_surv) <- c("est", "L95%", "U95%", "SE", "knots", "gamma0", "gamma1", "gamma2", "gamma3", "gamma4")
  return(summary_surv)
}

# This version includes the treatment effect and is more flexible
extract_surv_fit_updated <- function(surv_object) {
  
  # Combine core summary components into a tibble
  knots <- c(surv_object$knots, "")
  summary_surv <- cbind(surv_object$res.t, knots, surv_object$cov)
  summary_surv <- as_tibble(summary_surv)
  
  # --- 1️⃣ Identify spline (gamma) coefficients ---
  gamma_idx <- grep("^gamma", names(coef(surv_object)))
  n_gamma <- length(gamma_idx)
  gamma_names <- paste0("gamma", 0:(n_gamma - 1))
  
  # --- 2️⃣ Identify treatment (covariate) effects ---
  covar_names <- names(coef(surv_object))[-gamma_idx]
  n_covars <- length(covar_names)
  
  # --- 3️⃣ Build column names dynamically ---
  names(summary_surv) <- c("est", "L95%", "U95%", "SE", "knots",
                           gamma_names, covar_names)
  
  return(summary_surv)
}

# Assumes input is vector of mean, lower and upper bounds
format_continuous_result <- function(x, n_digits = 2) {
  formatted_result <- paste0(format(x[1], digits = n_digits, nsmall = n_digits), " (",
                             format(x[2], digits = n_digits, nsmall = n_digits), ", ",
                             format(x[3], digits = n_digits, nsmall = n_digits), ")")
  if(length(x) == 4) formatted_result <- paste0(formatted_result, " p=", 
                                                format(x[4], digits = n_digits, nsmall = n_digits))
  
  return(formatted_result)
}