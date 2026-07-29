############################################################
# EM analysis across simulation settings
############################################################

setwd("~/Box Sync/MelissaMartin_Dissertation/Paper3")
source("code/helper_functions.R")

############################################################
# 1. Q FUNCTION
############################################################

Q <- function(params, params_old, dat, n_tv=3) {
  n = nrow(dat)
  
  sens = params[1]
  spec = params[2]
  prev_arm1 = params[3]
  prev_arm2 = params[4]
  variance = params[5]
  a1_r_cor = params[6]
  a1_r_mis = params[7]
  a1_nr_cor_a2nr_1 = params[8]
  a1_nr_mis_a2nr_1 = params[9]
  a1_nr_cor_a2nr_2 = params[10]
  a1_nr_mis_a2nr_2 = params[11]
  b1_r_cor = params[12]
  b1_r_mis = params[13]
  b1_nr_cor_b2nr_1 = params[14]
  b1_nr_mis_b2nr_1 = params[15]
  b1_nr_cor_b2nr_2 = params[16]
  b1_nr_mis_b2nr_2 = params[17]
  
  sens_old = params_old[1]
  spec_old = params_old[2]
  prev_arm1_old = params_old[3]
  prev_arm2_old = params_old[4]
  variance_old = params_old[5]
  a1_r_cor_old = params_old[6]
  a1_r_mis_old = params_old[7]
  a1_nr_cor_a2nr_1_old = params_old[8]
  a1_nr_mis_a2nr_1_old = params_old[9]
  a1_nr_cor_a2nr_2_old = params_old[10]
  a1_nr_mis_a2nr_2_old = params_old[11]
  b1_r_cor_old = params_old[12]
  b1_r_mis_old = params_old[13]
  b1_nr_cor_b2nr_1_old = params_old[14]
  b1_nr_mis_b2nr_1_old = params_old[15]
  b1_nr_cor_b2nr_2_old = params_old[16]
  b1_nr_mis_b2nr_2_old = params_old[17]
  
  prev_old_vec = data.table::fifelse(dat$treatment_1==1, prev_arm1_old, prev_arm2_old)
  
  t = dat$n_tv_pos
  
  numerator = prev_old_vec * sens_old^t * (1-sens_old)^(n_tv - t)
  denominator = numerator + (1-prev_old_vec) * (1-spec_old)^t * spec_old^(n_tv - t)
  prob_condit_truepos = numerator / denominator
  
  idx_10  <- dat$treatment_1 == 1 & dat$treatment_2 == 0
  idx_11  <- dat$treatment_1 == 1 & dat$treatment_2 == 1
  idx_1m1 <- dat$treatment_1 == 1 & dat$treatment_2 == -1
  
  idx_20  <- dat$treatment_1 == 2 & dat$treatment_2 == 0
  idx_21  <- dat$treatment_1 == 2 & dat$treatment_2 == 1
  idx_2m1 <- dat$treatment_1 == 2 & dat$treatment_2 == -1
  
  mu_s1_old_vec <- numeric(nrow(dat))
  mu_s1_old_vec[idx_10]  <- a1_r_cor_old
  mu_s1_old_vec[idx_11]  <- a1_nr_mis_a2nr_1_old
  mu_s1_old_vec[idx_1m1] <- a1_nr_mis_a2nr_2_old
  mu_s1_old_vec[idx_20]  <- b1_r_cor_old
  mu_s1_old_vec[idx_21]  <- b1_nr_mis_b2nr_1_old
  mu_s1_old_vec[idx_2m1] <- b1_nr_mis_b2nr_2_old
  
  mu_s0_old_vec <- numeric(nrow(dat))
  mu_s0_old_vec[idx_10]  <- a1_r_mis_old
  mu_s0_old_vec[idx_11]  <- a1_nr_cor_a2nr_1_old
  mu_s0_old_vec[idx_1m1] <- a1_nr_cor_a2nr_2_old
  mu_s0_old_vec[idx_20]  <- b1_r_mis_old
  mu_s0_old_vec[idx_21]  <- b1_nr_cor_b2nr_1_old
  mu_s0_old_vec[idx_2m1] <- b1_nr_cor_b2nr_2_old
  
  dnorm_s1 = dnorm(dat$outcome, mean=mu_s1_old_vec, sd=sqrt(variance_old))
  dnorm_s0 = dnorm(dat$outcome, mean=mu_s0_old_vec, sd=sqrt(variance_old))
  
  w = (dnorm_s1 * prob_condit_truepos) /
    (dnorm_s1 * prob_condit_truepos + dnorm_s0 * (1-prob_condit_truepos))
  
  term1 = log(prev_arm1) * sum(w[which(dat$treatment_1==1)])
  term2 = log(1-prev_arm1) * sum((1 - w[which(dat$treatment_1==1)]))
  term3 = log(prev_arm2) * sum(w[which(dat$treatment_1==2)])
  term4 = log(1-prev_arm2) * sum((1 - w[which(dat$treatment_1==2)]))
  
  term5 = log(sens) * sum(w * t)
  term6 = log(1-sens) * sum(w * (n_tv - t))
  term7 = log(1-spec) * sum((1-w) * t)
  term8 = log(spec) * sum((1-w) * (n_tv - t))
  
  term9 = -(n/2) * log(variance)
  
  term10 = -(1/(2*variance)) * sum((dat$outcome[idx_10] - a1_r_cor)^2 * w[idx_10])
  term11 = -(1/(2*variance)) * sum((dat$outcome[idx_10] - a1_r_mis)^2 * (1 - w[idx_10]))
  
  term12 = -(1/(2*variance)) * sum((dat$outcome[idx_11] - a1_nr_cor_a2nr_1)^2 * (1 - w[idx_11]))
  term13 = -(1/(2*variance)) * sum((dat$outcome[idx_11] - a1_nr_mis_a2nr_1)^2 * w[idx_11])
  
  term14 = -(1/(2*variance)) * sum((dat$outcome[idx_1m1] - a1_nr_cor_a2nr_2)^2 * (1 - w[idx_1m1]))
  term15 = -(1/(2*variance)) * sum((dat$outcome[idx_1m1] - a1_nr_mis_a2nr_2)^2 * w[idx_1m1])
  
  term16 = -(1/(2*variance)) * sum((dat$outcome[idx_20] - b1_r_cor)^2 * w[idx_20])
  term17 = -(1/(2*variance)) * sum((dat$outcome[idx_20] - b1_r_mis)^2 * (1 - w[idx_20]))
  
  term18 = -(1/(2*variance)) * sum((dat$outcome[idx_21] - b1_nr_cor_b2nr_1)^2 * (1 - w[idx_21]))
  term19 = -(1/(2*variance)) * sum((dat$outcome[idx_21] - b1_nr_mis_b2nr_1)^2 * w[idx_21])
  
  term20 = -(1/(2*variance)) * sum((dat$outcome[idx_2m1] - b1_nr_cor_b2nr_2)^2 * (1 - w[idx_2m1]))
  term21 = -(1/(2*variance)) * sum((dat$outcome[idx_2m1] - b1_nr_mis_b2nr_2)^2 * w[idx_2m1])
  
  return(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 +
           term9 + term10 + term11 + term12 + term13 + term14 + term15 +
           term16 + term17 + term18 + term19 + term20 + term21)
}

############################################################
# 2. EM FUNCTION FOR ONE DATASET
############################################################

run_em <- function(dat_cur, n_tv = 3, tol = 1e-5, max_iter = 500) {
  
  params <- c(
    sensitivity = 0.7,
    specificity = 0.9,
    prevalence_arm1 = 0.7,
    prevalence_arm2 = 0.6,
    variance = 0.5,
    mu_1110 = 5,
    mu_1010 = 5,
    mu_1001 = 5,
    mu_1101 = 5,
    mu_1002 = 5,
    mu_1102 = 5,
    mu_2110 = 5,
    mu_2010 = 5,
    mu_2001 = 5,
    mu_2101 = 5,
    mu_2002 = 5,
    mu_2102 = 5
  )
  
  lower <- c(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, rep(-Inf, 12))
  upper <- c(1 - 1e-6, 1 - 1e-6, 1 - 1e-6, 1 - 1e-6, Inf, rep(Inf, 12))
  
  diff <- Inf
  iter <- 0
  converged <- FALSE
  last_value <- NA
  
  while (diff >= tol && iter < max_iter) {
    iter <- iter + 1
    params_old <- params
    
    o <- optim(
      par = params_old,
      fn = Q,
      params_old = params_old,
      dat = dat_cur,
      n_tv = n_tv,
      lower = lower,
      upper = upper,
      method = "L-BFGS-B",
      control = list(fnscale = -1)
    )
    
    params <- o$par
    diff <- max(abs(params - params_old))
    last_value <- o$value
  }
  
  if (diff < tol) converged <- TRUE
  
  return(list(
    est = params,
    converged = converged,
    iter = iter,
    diff = diff,
    Q_final = last_value
  ))
}

############################################################
# 3. SETTINGS
############################################################

n_sim <- 500
n_tv <- 3

base_path <- "~/Box Sync/MelissaMartin_Dissertation/Paper3"

settings <- expand.grid(
  prev_label = c("HIGH", "LOW"),
  sens_label = c("HIGH", "LOW"),
  spec_label = c("HIGH", "LOW"),
  EM_alloc = c(TRUE,FALSE)
)

all_recovery <- list()
all_dtr <- list()

############################################################
# 4. LOOP OVER SETTINGS
############################################################

for (s in 1:nrow(settings)) {
  
  prev_label <- settings$prev_label[s]
  sens_label <- settings$sens_label[s]
  spec_label <- settings$spec_label[s]
  EM_alloc <- settings$EM_alloc[s]
  
  path_use = ifelse(s %in% c(8,7,5,3,4),
                    paste0(base_path, "/nick_simulations/results"),
                    paste0(base_path, "/smart_simulations/results"))
  
  sim_path <- file.path(
    path_use,
    paste0("EM_alloc_", EM_alloc),
    paste0("prev", prev_label, "_sens", sens_label, "_spec", spec_label)
  )
  
  setting_name <- paste(
    ifelse(prev_label == "HIGH", "High prevalence", "Low prevalence"),
    ifelse(sens_label == "HIGH", "high sensitivity", "low sensitivity"),
    ifelse(spec_label == "HIGH", "high specificity", "low specificity"),
    sep = ", "
  )
  
  cat("\n====================================\n")
  cat(setting_name, "\n")
  cat("Path:", sim_path, "\n")
  cat("====================================\n")
  
  true_params <- c(
    sensitivity = ifelse(sens_label == "HIGH", 0.8, 0.5),
    specificity = ifelse(spec_label == "HIGH", 0.8, 0.5),
    prevalence_arm1 = ifelse(prev_label == "HIGH", 0.7, 0.4),
    prevalence_arm2 = ifelse(prev_label == "HIGH", 0.6, 0.3),
    variance = 1,
    mu_1110 = 10,
    mu_1010 = 7,
    mu_1001 = 9,
    mu_1101 = 8,
    mu_1002 = 8,
    mu_1102 = 7,
    mu_2110 = 8,
    mu_2010 = 6,
    mu_2001 = 8,
    mu_2101 = 7,
    mu_2002 = 7,
    mu_2102 = 6
  )
  
  results <- vector("list", n_sim)
  
  for (b in 1:n_sim) {
    if (b %% 50 == 0) cat("Simulation", b, "\n")
    
    dat_file <- ifelse(s==5,
                       file.path(sim_path, paste0("sim", b+500), "dat_upd.RData"),
                       file.path(sim_path, paste0("sim", b), "dat_upd.RData"))
    
    results[[b]] <- tryCatch({
      load(dat_file)  # loads dat
      out <- run_em(dat_cur = dat, n_tv = n_tv)
      names(out$est) <- names(true_params)
      out
    }, error = function(e) {
      message("Error in setting ", s, ", sim ", b, ": ", e$message)
      list(
        est = setNames(rep(NA, length(true_params)), names(true_params)),
        converged = FALSE,
        iter = NA,
        diff = NA,
        Q_final = NA
      )
    })
  }
  
  save(results, file = file.path(sim_path, "results_EMalg_analysis.RData"))
  
  ############################################################
  # 5. ESTIMATES DATA FRAME
  ############################################################
  
  est_mat <- do.call(rbind, lapply(results, function(x) x$est))
  est_df <- as.data.frame(est_mat)
  
  conv_vec <- sapply(results, function(x) x$converged)
  iter_vec <- sapply(results, function(x) x$iter)
  
  cat("Convergence rate:", mean(conv_vec, na.rm = TRUE), "\n")
  print(summary(iter_vec))
  
  est_conv <- est_df[conv_vec, , drop = FALSE]
  
  ############################################################
  # 6. PARAMETER RECOVERY
  ############################################################
  
  recovery_summary <- data.frame(
    simulation_setting = setting_name,
    parameter = names(true_params),
    rel_bias = NA,
    rmse = NA,
    convergence_rate = mean(conv_vec, na.rm = TRUE)
  )
  
  for (j in seq_along(true_params)) {
    p <- names(true_params)[j]
    est_j <- est_conv[[p]]
    true_j <- true_params[p]
    
    recovery_summary$rel_bias[j] <- mean((est_j - true_j) / true_j, na.rm = TRUE)
    recovery_summary$rmse[j] <- sqrt(mean((est_j - true_j)^2, na.rm = TRUE))
  }
  
  all_recovery[[s]] <- recovery_summary
  
  ############################################################
  # 7. TRUE IDEAL DTR MEANS
  ############################################################
  
  true_dtr <- c(
    dtr_1_1  = true_params["prevalence_arm1"] * true_params["mu_1110"] +
      (1 - true_params["prevalence_arm1"]) * true_params["mu_1001"],
    
    dtr_1_m1 = true_params["prevalence_arm1"] * true_params["mu_1110"] +
      (1 - true_params["prevalence_arm1"]) * true_params["mu_1002"],
    
    dtr_2_1  = true_params["prevalence_arm2"] * true_params["mu_2110"] +
      (1 - true_params["prevalence_arm2"]) * true_params["mu_2001"],
    
    dtr_2_m1 = true_params["prevalence_arm2"] * true_params["mu_2110"] +
      (1 - true_params["prevalence_arm2"]) * true_params["mu_2002"]
  )
  
  true_dtr <- as.numeric(true_dtr)
  names(true_dtr) <- c("dtr_1_1", "dtr_1_m1", "dtr_2_1", "dtr_2_m1")
  
  ############################################################
  # 8. ESTIMATED IDEAL DTR MEANS
  ############################################################
  
  dtr_est <- data.frame(
    dtr_1_1 = est_conv$prevalence_arm1 * est_conv$mu_1110 +
      (1 - est_conv$prevalence_arm1) * est_conv$mu_1001,
    
    dtr_1_m1 = est_conv$prevalence_arm1 * est_conv$mu_1110 +
      (1 - est_conv$prevalence_arm1) * est_conv$mu_1002,
    
    dtr_2_1 = est_conv$prevalence_arm2 * est_conv$mu_2110 +
      (1 - est_conv$prevalence_arm2) * est_conv$mu_2001,
    
    dtr_2_m1 = est_conv$prevalence_arm2 * est_conv$mu_2110 +
      (1 - est_conv$prevalence_arm2) * est_conv$mu_2002
  )
  
  dtr_summary <- data.frame(
    simulation_setting = setting_name,
    DTR = names(true_dtr),
    rel_bias = NA,
    rmse = NA,
    convergence_rate = mean(conv_vec, na.rm = TRUE)
  )
  
  for (j in seq_along(true_dtr)) {
    d <- names(true_dtr)[j]
    est_j <- dtr_est[[d]]
    true_j <- true_dtr[d]
    
    dtr_summary$rel_bias[j] <- mean((est_j - true_j) / true_j, na.rm = TRUE)
    dtr_summary$rmse[j] <- sqrt(mean((est_j - true_j)^2, na.rm = TRUE))
  }
  
  all_dtr[[s]] <- dtr_summary
}

############################################################
# 9. COMBINE AND SAVE
############################################################

recovery_all <- do.call(rbind, all_recovery)
dtr_all <- do.call(rbind, all_dtr)

save(
  recovery_all,
  dtr_all,
  file = file.path(base_path, "smart_simulations/recovery_all_settings.RData")
)

write.csv(
  recovery_all,
  file = file.path(base_path, "smart_simulations/parameter_recovery_all_settings.csv"),
  row.names = FALSE
)

write.csv(
  dtr_all,
  file = file.path(base_path, "smart_simulations/dtr_recovery_all_settings.csv"),
  row.names = FALSE
)