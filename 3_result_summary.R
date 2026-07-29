setwd("~/Box Sync/MelissaMartin_Dissertation/Paper3")
source("code/helper_functions.R")
library(dplyr)
library(tidyr)


############################################################
# 5. ESTIMATES DATA FRAME
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
all_plot_params <- list()
all_plot_dtr <- list()

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
  
  load(file.path(sim_path, "results_EMalg_analysis.RData"))
  
  est_mat <- do.call(rbind, lapply(results, function(x) x$est))
  est_df <- as.data.frame(est_mat)
  
  conv_vec <- sapply(results, function(x) x$converged)
  iter_vec <- sapply(results, function(x) x$iter)
  
  cat("Convergence rate:", mean(conv_vec, na.rm = TRUE), "\n")
  print(summary(iter_vec))
  
  #est_conv <- est_df[conv_vec, , drop = FALSE]
  
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
    est_j <- est_df[[p]]
    true_j <- true_params[p]
    
    rel_bias_vec <- (est_j - true_j) / true_j
    error_vec <- est_j - true_j
    
    recovery_summary$rel_bias[j] <- mean(rel_bias_vec, na.rm = TRUE) * 100
    recovery_summary$rmse[j] <- sqrt(mean(error_vec^2, na.rm = TRUE))
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
    dtr_1_1 = est_df$prevalence_arm1 * est_df$mu_1110 +
      (1 - est_df$prevalence_arm1) * est_df$mu_1001,
    
    dtr_1_m1 = est_df$prevalence_arm1 * est_df$mu_1110 +
      (1 - est_df$prevalence_arm1) * est_df$mu_1002,
    
    dtr_2_1 = est_df$prevalence_arm2 * est_df$mu_2110 +
      (1 - est_df$prevalence_arm2) * est_df$mu_2001,
    
    dtr_2_m1 = est_df$prevalence_arm2 * est_df$mu_2110 +
      (1 - est_df$prevalence_arm2) * est_df$mu_2002
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
    
    rel_bias_vec <- (est_j - true_j) / true_j
    error_vec <- est_j - true_j
    
    dtr_summary$rel_bias[j] <- mean(rel_bias_vec, na.rm = TRUE) * 100
    dtr_summary$rmse[j] <- sqrt(mean((est_j - true_j)^2, na.rm = TRUE))
  }
  
  all_dtr[[s]] <- dtr_summary
  
  if (prev_label == "HIGH" & EM_alloc == TRUE &
      ((sens_label == "HIGH" & spec_label == "HIGH") |
       (sens_label == "LOW"  & spec_label == "LOW"))) {
    
    setting_short <- ifelse(
      sens_label == "HIGH" & spec_label == "HIGH",
      "HHH",
      "HLL"
    )
    
    est_conv <- est_df
    est_conv$sim <- 1:nrow(est_conv)
    
    est_long <- pivot_longer(
      est_conv,
      cols = -sim,
      names_to = "parameter",
      values_to = "estimate"
    )
    
    est_long$true <- true_params[est_long$parameter]
    est_long$rel_bias <- (est_long$estimate - est_long$true) / est_long$true
    
    assign(paste0("est_long_", setting_short), est_long)
    
    dtr_conv <- dtr_est
    dtr_conv$sim <- 1:nrow(dtr_conv)
    
    dtr_long <- pivot_longer(
      dtr_conv,
      cols = -sim,
      names_to = "DTR",
      values_to = "estimate"
    )
    
    dtr_long$true <- true_dtr[dtr_long$DTR]
    dtr_long$rel_bias <- (dtr_long$estimate - dtr_long$true) / dtr_long$true
    
    assign(paste0("dtr_long_", setting_short), dtr_long)
  }
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

param_plot_all <- rbindlist(all_plot_params)
dtr_plot_all <- rbindlist(all_plot_dtr)
recovery_summary_all = rbindlist(all_recovery)


### figures for the HHH and HLL settings
library(dplyr)
library(ggplot2)

param_labels <- c(
  sensitivity = expression(Se),
  specificity = expression(Sp),
  prevalence_arm1 = expression(pi[1]),
  prevalence_arm2 = expression(pi[-1]),
  variance = expression(sigma^2),
  
  # A1 = 1
  mu_1110 = expression(mu[paste(1, ",", 1)]^{paste(1, ",", 0)}),
  mu_1001 = expression(mu[paste(0, ",", 0)]^{paste(1, ",", 1)}),
  mu_1002 = expression(mu[paste(0, ",", 0)]^{paste(1, ",", -1)}),
  mu_1010 = expression(mu[paste(0, ",", 1)]^{paste(1, ",", 0)}),
  mu_1101 = expression(mu[paste(1, ",", 0)]^{paste(1, ",", 1)}),
  mu_1102 = expression(mu[paste(1, ",", 0)]^{paste(1, ",", -1)}),
  
  # A1 = -1
  mu_2110 = expression(mu[paste(1, ",", 1)]^{paste(-1, ",", 0)}),
  mu_2001 = expression(mu[paste(0, ",", 0)]^{paste(-1, ",", 1)}),
  mu_2002 = expression(mu[paste(0, ",", 0)]^{paste(-1, ",", -1)}),
  mu_2010 = expression(mu[paste(0, ",", 1)]^{paste(-1, ",", 0)}),
  mu_2101 = expression(mu[paste(1, ",", 0)]^{paste(-1, ",", 1)}),
  mu_2102 = expression(mu[paste(1, ",", 0)]^{paste(-1, ",", -1)})
)

param_order <- c(
  # A1 = 1
  "mu_1110",
  "mu_1001",
  "mu_1002",
  "mu_1010",
  "mu_1101",
  "mu_1102",
  
  # A1 = -1
  "mu_2110",
  "mu_2001",
  "mu_2002",
  "mu_2010",
  "mu_2101",
  "mu_2102",
  
  # other parameters
  "variance",
  "sensitivity",
  "specificity",
  "prevalence_arm1",
  "prevalence_arm2"
)

est_long_HHH$parameter <- factor(est_long_HHH$parameter, levels = param_order)
est_long_HLL$parameter <- factor(est_long_HLL$parameter, levels = param_order)

p_param_HHH = ggplot(est_long_HHH, aes(x = parameter, y = rel_bias)) +
  geom_boxplot(fill = "gray90") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = param_labels) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  ) +
  labs(
    x = "Parameter",
    y = "Relative bias"
  )

ggsave("~/Box Sync/MelissaMartin_Dissertation/Paper3/figures/parameter_relbias_HHH.png", p_param_HHH, width=12, height=8, dpi=300)

p_param_HLL = ggplot(est_long_HLL, aes(x = parameter, y = rel_bias)) +
  geom_boxplot(fill = "gray90") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = param_labels) +
  theme_bw() +
  theme(
    #axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
  ) +
  labs(
    x = "Parameter",
    y = "Relative bias",
  )

ggsave("~/Box Sync/MelissaMartin_Dissertation/Paper3/figures/parameter_relbias_HLL.png", p_param_HLL, width=12, height=8, dpi=300)


dtr_labels <- c(
  dtr_1_1  = "(1, 1)",
  dtr_1_m1 = "(1, -1)",
  dtr_2_1  = "(-1, 1)",
  dtr_2_m1 = "(-1, -1)"
)

dtr_order <- names(dtr_labels)

dtr_long_HHH$DTR <- factor(dtr_long_HHH$DTR, levels = dtr_order)
dtr_long_HLL$DTR <- factor(dtr_long_HLL$DTR, levels = dtr_order)

p_dtr_HHH = ggplot(dtr_long_HHH, aes(x = DTR, y = rel_bias)) +
  geom_boxplot(fill = "gray90") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = dtr_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
    ) +
  labs(
    x = "Ideal DTR",
    y = "Relative bias"
  )

ggsave("~/Box Sync/MelissaMartin_Dissertation/Paper3/figures/dtr_relbias_HHH.png", p_dtr_HHH, width=8.5, height=5, dpi=300)


p_dtr_HLL = ggplot(dtr_long_HLL, aes(x = DTR, y = rel_bias)) +
  geom_boxplot(fill = "gray90") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_discrete(labels = dtr_labels) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15)
    ) +
  labs(
    x = "Ideal DTR",
    y = "Relative bias"
  )

ggsave("~/Box Sync/MelissaMartin_Dissertation/Paper3/figures/dtr_relbias_HLL.png", p_dtr_HLL, width=8.5, height=5, dpi=300)


